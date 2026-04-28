"""fetch-references — idempotent downloader for Taiji-agent reference data.

Reads `reference_manifest.yml` (sibling file) for per-genome file specs,
downloads anything missing into <output_dir>/, gunzips on the fly when the
manifest says so, and (optionally) builds a samtools faidx index.

Idempotency rules:
- A target file present at the expected path with size within +/-5% of the
  manifest's approx_size_mb is treated as already-downloaded and skipped.
- --force re-downloads everything for the chosen genome.
- The skill never deletes files it didn't write itself (no destructive ops).

Library API
-----------

    from fetch import fetch_genome, check_genome
    result = fetch_genome("hg38", output_dir="./dependencies_data")
    # result.files: {kind: Path}     e.g. {"fasta": "./dependencies_data/hg38/genome.fa", ...}

CLI
---

    python fetch.py --list
    python fetch.py --genome hg38 --output dependencies_data/
    python fetch.py --genome hg38 --output dependencies_data/ --check
    python fetch.py --genome hg38 --output dependencies_data/ --force
    python fetch.py --genome hg38 --output dependencies_data/ --update-genomes-yml
    python fetch.py --genome hg38 --output dependencies_data/ --dry-run
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

# Motif sources supported by the manifest. New ones added as `motif_<key>:`
# entries in reference_manifest.yml. Default is set by --motif-source.
MOTIF_SOURCES = ("cisbp", "hocomoco")
DEFAULT_MOTIF_SOURCE = "cisbp"

try:
    import yaml
except ImportError:
    print("fetch-references: pyyaml is required. "
          "`pip install pyyaml` or `micromamba install pyyaml`.", file=sys.stderr)
    sys.exit(2)


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_MANIFEST = SCRIPT_DIR / "reference_manifest.yml"


# ---------------------------------------------------------------------------
# Soft attach to workflow-log
# ---------------------------------------------------------------------------


def _attach_log():
    """Best-effort attach to workflow-log skill's active run."""
    try:
        sibling = SCRIPT_DIR.parent.parent / "workflow-log" / "scripts"
        if sibling.exists() and str(sibling) not in sys.path:
            sys.path.insert(0, str(sibling))
        from log import TaijiLog  # type: ignore[import-not-found]
        return TaijiLog.attach_active(work_dir=Path.cwd())
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass
class FileResult:
    kind: str            # "fasta" | "gtf" | "motif"
    target: Path
    status: str          # "present" | "downloaded" | "skipped" | "failed" | "dry-run"
    bytes_on_disk: int = 0
    error: str | None = None


@dataclass
class FetchResult:
    genome: str
    output_dir: Path
    files: dict[str, FileResult] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return all(f.status in ("present", "downloaded", "dry-run")
                   for f in self.files.values())

    def to_dict(self) -> dict:
        return {
            "genome": self.genome,
            "output_dir": str(self.output_dir),
            "ok": self.ok,
            "files": {k: {"target": str(v.target), "status": v.status,
                          "bytes_on_disk": v.bytes_on_disk, "error": v.error}
                      for k, v in self.files.items()},
        }


# ---------------------------------------------------------------------------
# Manifest loading
# ---------------------------------------------------------------------------


def load_manifest(path: Path = DEFAULT_MANIFEST) -> dict:
    if not path.exists():
        raise SystemExit(f"manifest not found: {path}")
    with path.open() as fh:
        return yaml.safe_load(fh) or {}


def list_genomes(manifest: dict) -> None:
    print("Genomes available in reference_manifest.yml:\n")
    for name, entry in manifest.items():
        print(f"  {name:6s}  {entry.get('description', '')}")
        for kind, spec in (entry.get("files") or {}).items():
            sz = spec.get("approx_size_mb", "?")
            print(f"           - {kind:6s} ~{sz:>5} MB  {spec['url']}")
        print()


# ---------------------------------------------------------------------------
# Download primitives
# ---------------------------------------------------------------------------


def _within_tolerance(actual_bytes: int, expected_mb: float | None,
                      pct: float = 5.0) -> bool:
    """Check actual file size is within +/-pct% of expected_mb."""
    if expected_mb is None:
        return True   # no expected size declared -> trust whatever's there
    expected = expected_mb * 1024 * 1024
    lo = expected * (1 - pct / 100)
    hi = expected * (1 + pct / 100)
    return lo <= actual_bytes <= hi


def _human_bytes(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n = n / 1024.0
    return f"{n:.1f} TB"


def _download(url: str, dest: Path, gunzip: bool, log_prefix: str = "") -> int:
    """Stream `url` to `dest`. If gunzip=True, transparently inflate the
    download stream into `dest`. Returns bytes written to disk.

    Tries wget -> curl -> urllib in that order. wget/curl honor the system
    CA trust store, which matters on cluster nodes where the conda env's
    OpenSSL bundle may be missing/stale and Python alone can't validate
    HTTPS. Override the order with TAIJI_FETCH_DOWNLOADER=wget|curl|urllib.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"{log_prefix}downloading {url}", file=sys.stderr)
    t0 = time.time()
    try:
        bytes_written = _do_download(url, tmp, gunzip, log_prefix, t0)
        tmp.rename(dest)
    except Exception:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass
        raise
    elapsed = time.time() - t0
    rate = bytes_written / (elapsed + 1e-9)
    print(f"{log_prefix}  -> {dest}  {_human_bytes(bytes_written)} "
          f"in {elapsed:.1f}s ({_human_bytes(int(rate))}/s)", file=sys.stderr)
    return bytes_written


def _do_download(url: str, tmp: Path, gunzip: bool, log_prefix: str, t0: float) -> int:
    override = os.environ.get("TAIJI_FETCH_DOWNLOADER", "").lower()
    if override in ("wget", "curl", "urllib"):
        tools = [override]
    else:
        tools = []
        if shutil.which("wget"):
            tools.append("wget")
        if shutil.which("curl"):
            tools.append("curl")
        tools.append("urllib")

    last_err: Exception | None = None
    for tool in tools:
        try:
            print(f"{log_prefix}  using {tool}", file=sys.stderr)
            if tool == "wget":
                return _via_subprocess(["wget", "-q", "-O", "-", url], tmp, gunzip)
            if tool == "curl":
                return _via_subprocess(["curl", "-fsSL", url], tmp, gunzip)
            return _via_urllib(url, tmp, gunzip, log_prefix, t0)
        except Exception as e:
            last_err = e
            print(f"{log_prefix}  {tool} failed: {e}", file=sys.stderr)
    assert last_err is not None
    raise last_err


def _via_subprocess(fetch_cmd: list[str], tmp: Path, gunzip: bool) -> int:
    """Run wget/curl, optionally piping through gunzip into tmp. Bytes written = tmp size."""
    with tmp.open("wb") as out:
        if gunzip:
            p1 = subprocess.Popen(fetch_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            assert p1.stdout is not None
            p2 = subprocess.Popen(["gunzip", "-c"], stdin=p1.stdout, stdout=out,
                                  stderr=subprocess.PIPE)
            p1.stdout.close()
            _, err2 = p2.communicate()
            _, err1 = p1.communicate()
            if p1.returncode != 0:
                raise RuntimeError(
                    f"{fetch_cmd[0]} exit {p1.returncode}: {err1.decode(errors='replace')[:300]}")
            if p2.returncode != 0:
                raise RuntimeError(
                    f"gunzip exit {p2.returncode}: {err2.decode(errors='replace')[:300]}")
        else:
            r = subprocess.run(fetch_cmd, stdout=out, stderr=subprocess.PIPE)
            if r.returncode != 0:
                raise RuntimeError(
                    f"{fetch_cmd[0]} exit {r.returncode}: {r.stderr.decode(errors='replace')[:300]}")
    return tmp.stat().st_size


def _via_urllib(url: str, tmp: Path, gunzip: bool, log_prefix: str, t0: float) -> int:
    bytes_written = 0
    last_report = t0
    with urllib.request.urlopen(url, timeout=60) as resp:
        src = gzip.GzipFile(fileobj=resp) if gunzip else resp  # type: ignore[arg-type]
        with tmp.open("wb") as out:
            while True:
                chunk = src.read(1 << 20)   # 1 MB chunks
                if not chunk:
                    break
                out.write(chunk)
                bytes_written += len(chunk)
                now = time.time()
                if now - last_report > 5:
                    rate = bytes_written / (now - t0 + 1e-9)
                    print(f"{log_prefix}  ... {_human_bytes(bytes_written)} "
                          f"({_human_bytes(int(rate))}/s)", file=sys.stderr)
                    last_report = now
    return bytes_written


def _build_fai(fasta: Path) -> bool:
    """Run `samtools faidx` if samtools is on PATH. Returns True on success."""
    if shutil.which("samtools") is None:
        print(f"  NOTE: samtools not on PATH; skipping faidx for {fasta.name}. "
              "Run `samtools faidx <path>` manually after install.",
              file=sys.stderr)
        return False
    try:
        subprocess.run(["samtools", "faidx", str(fasta)],
                       check=True, capture_output=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"  WARNING: samtools faidx failed: {e.stderr.decode()[:200]}",
              file=sys.stderr)
        return False


# ---------------------------------------------------------------------------
# Motif-source dispatch
# ---------------------------------------------------------------------------


def _resolve_motif_source(files_spec: dict, motif_source: str) -> dict:
    """Filter the per-genome files dict so only the chosen motif source
    remains, renamed to the bare `motif` key for downstream reporting.

    The manifest carries parallel `motif_cisbp:` and `motif_hocomoco:` blocks
    per genome; this picks one and discards the other so the rest of the
    pipeline doesn't have to know about source-specific keys.
    """
    out: dict = {}
    chosen_key = f"motif_{motif_source}"
    for k, v in files_spec.items():
        if not k.startswith("motif_"):
            out[k] = v
            continue
        if k == chosen_key:
            out["motif"] = v
    if "motif" not in out:
        avail = sorted(k.removeprefix("motif_") for k in files_spec
                       if k.startswith("motif_"))
        raise SystemExit(
            f"motif source '{motif_source}' not declared in manifest for "
            f"this genome. Available: {avail}. Add a `motif_{motif_source}:` "
            f"block to reference_manifest.yml or pick a different source."
        )
    return out


# ---------------------------------------------------------------------------
# Tarball extraction (for sources like CIS-BP that ship inside MEME Suite's
# ~30 MB motif-database tarball — we want one .meme file out of it)
# ---------------------------------------------------------------------------


def _extract_from_tar(tar_path: Path, member_path: str, dest: Path,
                      log_prefix: str = "") -> int:
    """Extract a single member from a (compressed) tar file to `dest`.

    Atomic via .part rename. Returns bytes written. Raises if member missing.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"{log_prefix}extracting {member_path} from {tar_path.name}",
          file=sys.stderr)
    with tarfile.open(tar_path, "r:*") as tf:
        try:
            member = tf.getmember(member_path)
        except KeyError:
            available = [m.name for m in tf.getmembers()
                         if member_path.split("/")[-1] in m.name][:10]
            raise RuntimeError(
                f"member '{member_path}' not in {tar_path.name}. "
                f"Closest matches: {available}"
            )
        src = tf.extractfile(member)
        if src is None:
            raise RuntimeError(f"member '{member_path}' is not a regular file")
        with tmp.open("wb") as out:
            shutil.copyfileobj(src, out, length=1 << 20)
    tmp.rename(dest)
    return dest.stat().st_size


def _fetch_to_cache(url: str, cache_path: Path, log_prefix: str = "") -> Path:
    """Download `url` to `cache_path` if not already there. No expected-size
    check — caller's manifest entry vouches for it."""
    if cache_path.exists() and cache_path.stat().st_size > 0:
        print(f"{log_prefix}cache hit: {cache_path}", file=sys.stderr)
        return cache_path
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    # _download() already does atomic .part rename + wget/curl/urllib fallback.
    _download(url, cache_path, gunzip=False, log_prefix=log_prefix)
    return cache_path


# ---------------------------------------------------------------------------
# Stable motif symlink
# ---------------------------------------------------------------------------


def _refresh_motifs_symlink(motif_target: Path) -> Path | None:
    """Create or update <genome>/motifs.meme as a relative symlink to the
    just-fetched motif file. Lets downstream Taiji templates reference a
    stable filename regardless of which motif source is in play.

    Returns the symlink path on success, None if the target's parent dir
    doesn't exist (shouldn't happen post-fetch).
    """
    if not motif_target.exists():
        return None
    link = motif_target.parent / "motifs.meme"
    # Replace if present; preserve absent.
    if link.is_symlink() or link.exists():
        try:
            link.unlink()
        except OSError:
            pass
    try:
        link.symlink_to(motif_target.name)   # relative target name
    except OSError as e:
        print(f"  WARN: could not create motifs.meme symlink: {e}",
              file=sys.stderr)
        return None
    return link


# ---------------------------------------------------------------------------
# Per-file orchestration
# ---------------------------------------------------------------------------


def _process_file(kind: str, spec: dict, output_dir: Path, *,
                  force: bool, dry_run: bool, check_only: bool) -> FileResult:
    target = output_dir / spec["target"]
    expected_mb = spec.get("approx_size_mb")
    gunzip = bool(spec.get("gunzip", False))
    url = spec["url"]

    # Check if already present and the right size.
    if target.exists() and not force:
        size = target.stat().st_size
        if _within_tolerance(size, expected_mb):
            return FileResult(kind=kind, target=target,
                              status="present", bytes_on_disk=size)
        if check_only:
            return FileResult(kind=kind, target=target,
                              status="present",
                              bytes_on_disk=size,
                              error=f"size {_human_bytes(size)} outside "
                                    f"+/-5% of expected {expected_mb} MB")
        # else: fall through to re-download

    if check_only:
        return FileResult(kind=kind, target=target,
                          status="skipped", bytes_on_disk=0,
                          error="not present (use without --check to download)")

    if dry_run:
        return FileResult(kind=kind, target=target,
                          status="dry-run", bytes_on_disk=0,
                          error=f"would download {url}")

    # Tarball-extraction path: download URL once into <output>/_cache/, then
    # extract just the named member into the target.
    extract_path = spec.get("extract_from_tar")
    try:
        if extract_path:
            cache_dir = output_dir / "_cache"
            cache_path = cache_dir / Path(url).name
            _fetch_to_cache(url, cache_path, log_prefix=f"  [{kind}] ")
            bytes_written = _extract_from_tar(cache_path, extract_path, target,
                                              log_prefix=f"  [{kind}] ")
        else:
            bytes_written = _download(url, target, gunzip=gunzip,
                                      log_prefix=f"  [{kind}] ")
    except Exception as e:
        return FileResult(kind=kind, target=target,
                          status="failed", error=str(e)[:240])

    # Optional: build samtools faidx.
    if spec.get("build_fai"):
        _build_fai(target)

    return FileResult(kind=kind, target=target,
                      status="downloaded", bytes_on_disk=bytes_written)


# ---------------------------------------------------------------------------
# Top-level fetch entry point (library API)
# ---------------------------------------------------------------------------


def fetch_genome(genome: str, output_dir: str | Path, *,
                 manifest: dict | None = None,
                 motif_source: str = DEFAULT_MOTIF_SOURCE,
                 force: bool = False, dry_run: bool = False,
                 check_only: bool = False) -> FetchResult:
    if manifest is None:
        manifest = load_manifest()
    if genome not in manifest:
        raise SystemExit(
            f"unknown genome '{genome}'. Available: {sorted(manifest)}")
    if motif_source not in MOTIF_SOURCES:
        raise SystemExit(
            f"unknown --motif-source '{motif_source}'. "
            f"Available: {list(MOTIF_SOURCES)}")
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    result = FetchResult(genome=genome, output_dir=output_dir)
    spec = manifest[genome]
    print(f"[fetch-references] genome={genome}  output_dir={output_dir}",
          file=sys.stderr)
    print(f"[fetch-references] motif source: {motif_source}", file=sys.stderr)
    print(f"[fetch-references] description: {spec.get('description', '')}",
          file=sys.stderr)
    files_spec = _resolve_motif_source(spec.get("files") or {}, motif_source)
    for kind, fspec in files_spec.items():
        result.files[kind] = _process_file(
            kind, fspec, output_dir,
            force=force, dry_run=dry_run, check_only=check_only,
        )
    # Refresh the stable motifs.meme symlink so Taiji templates referencing
    # <genome>/motifs.meme don't have to know which source was picked.
    motif_fr = result.files.get("motif")
    if motif_fr and motif_fr.status in ("downloaded", "present") and not dry_run:
        _refresh_motifs_symlink(motif_fr.target)
    return result


def check_genome(genome: str, output_dir: str | Path,
                 manifest: dict | None = None) -> FetchResult:
    """Convenience wrapper: --check semantics."""
    return fetch_genome(genome, output_dir, manifest=manifest, check_only=True)


# ---------------------------------------------------------------------------
# Optional: write the resolved paths into build-taiji-input's genomes.yml
# ---------------------------------------------------------------------------


def update_genomes_yml(result: FetchResult, genomes_yml: Path) -> None:
    """Edit build-taiji-input/assets/genomes.yml in place to point at the
    files we just downloaded. Only updates the entry for `result.genome`.
    Other entries are preserved verbatim."""
    if not genomes_yml.exists():
        print(f"  WARN: {genomes_yml} not found; skipping --update-genomes-yml",
              file=sys.stderr)
        return
    with genomes_yml.open() as fh:
        existing = yaml.safe_load(fh) or {}
    fasta = result.files.get("fasta")
    gtf   = result.files.get("gtf")
    if not fasta or not gtf:
        print("  WARN: fasta/gtf missing in fetch result; "
              "not updating genomes.yml", file=sys.stderr)
        return
    existing.setdefault(result.genome, {})
    existing[result.genome]["fasta"] = str(fasta.target.resolve())
    existing[result.genome]["gtf"]   = str(gtf.target.resolve())
    # Preserve any prior `hic` entry.
    with genomes_yml.open("w") as fh:
        yaml.safe_dump(existing, fh, sort_keys=False)
    print(f"[fetch-references] updated {genomes_yml} "
          f"with new {result.genome} paths", file=sys.stderr)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _print_report(result: FetchResult) -> None:
    print()
    print(f"=== fetch-references report — {result.genome} ===")
    width = max(len(k) for k in result.files) if result.files else 6
    for kind, fr in result.files.items():
        size = _human_bytes(fr.bytes_on_disk) if fr.bytes_on_disk else "—"
        line = f"  {kind:<{width}}  [{fr.status:<10}]  {fr.target}  ({size})"
        print(line)
        if fr.error:
            print(f"  {' ' * width}  -> {fr.error}")
    print(f"\noverall: {'OK' if result.ok else 'INCOMPLETE'}")


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST,
                   help="Path to reference_manifest.yml. Default: sibling file.")
    p.add_argument("--genome", help="Genome key (hg38, hg19, mm10, mm39, ...).")
    p.add_argument("--output", type=Path, default=Path("dependencies_data"),
                   help="Output directory. Default: dependencies_data/")
    p.add_argument("--list", action="store_true",
                   help="List available genomes and exit.")
    p.add_argument("--check", action="store_true",
                   help="Report what's present/missing; don't download.")
    p.add_argument("--force", action="store_true",
                   help="Re-download even if files already exist.")
    p.add_argument("--dry-run", action="store_true",
                   help="Resolve URLs and targets; don't download.")
    p.add_argument("--update-genomes-yml", action="store_true",
                   help="After download, edit "
                        "skills/build-taiji-input/assets/genomes.yml in place "
                        "with the resolved fasta/gtf paths.")
    p.add_argument("--motif-source", choices=list(MOTIF_SOURCES),
                   default=DEFAULT_MOTIF_SOURCE,
                   help=f"Which motif database to fetch. Default: "
                        f"{DEFAULT_MOTIF_SOURCE} (CIS-BP 2.00 — JASPAR + "
                        "TRANSFAC + ENCODE + SELEX-seq aggregated, keyed by "
                        "TF gene symbol). Alternative: hocomoco (HOCOMOCO v11).")
    p.add_argument("--json", action="store_true",
                   help="Print machine-readable JSON report instead of text.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    manifest = load_manifest(args.manifest)
    if args.list:
        list_genomes(manifest)
        return 0
    if not args.genome:
        print("--genome is required (or pass --list to see options)",
              file=sys.stderr)
        return 2

    result = fetch_genome(args.genome, args.output, manifest=manifest,
                          motif_source=args.motif_source,
                          force=args.force, dry_run=args.dry_run,
                          check_only=args.check)

    if args.json:
        print(json.dumps(result.to_dict(), indent=2))
    else:
        _print_report(result)

    # Best-effort: update build-taiji-input/assets/genomes.yml if asked.
    if args.update_genomes_yml and result.ok and not args.dry_run and not args.check:
        update_genomes_yml(
            result,
            Path(__file__).resolve().parent.parent.parent
            / "build-taiji-input" / "assets" / "genomes.yml",
        )

    # Best-effort: log to workflow-log if active.
    log = _attach_log()
    if log is not None:
        try:
            log.append_custom(
                "fetch_references",
                f"Stage 0 — fetch references ({result.genome})",
                result.to_dict(),
                status="pass" if result.ok else "fail",
            )
        except Exception:
            pass

    return 0 if result.ok else 1


if __name__ == "__main__":
    sys.exit(main())
