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
import time
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

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
    download stream into `dest`. Returns bytes written to disk."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    bytes_in_stream = 0
    bytes_written = 0
    t0 = time.time()
    last_report = t0
    print(f"{log_prefix}downloading {url}", file=sys.stderr)
    try:
        with urllib.request.urlopen(url, timeout=60) as resp:
            if gunzip:
                # Stream-decompress: wrap the response in a GzipFile.
                src = gzip.GzipFile(fileobj=resp)  # type: ignore[arg-type]
            else:
                src = resp
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

    # Download.
    try:
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
                 force: bool = False, dry_run: bool = False,
                 check_only: bool = False) -> FetchResult:
    if manifest is None:
        manifest = load_manifest()
    if genome not in manifest:
        raise SystemExit(
            f"unknown genome '{genome}'. Available: {sorted(manifest)}")
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    result = FetchResult(genome=genome, output_dir=output_dir)
    spec = manifest[genome]
    print(f"[fetch-references] genome={genome}  output_dir={output_dir}",
          file=sys.stderr)
    print(f"[fetch-references] description: {spec.get('description', '')}",
          file=sys.stderr)
    for kind, fspec in (spec.get("files") or {}).items():
        result.files[kind] = _process_file(
            kind, fspec, output_dir,
            force=force, dry_run=dry_run, check_only=check_only,
        )
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
