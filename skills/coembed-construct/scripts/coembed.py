"""Orchestrator for coembed-construct.

Thin Python CLI that:
  1. Validates the two input objects exist (one RNA, one ATAC .rds).
  2. Optionally gates on detect-dataset-type for each input, mostly for
     warning-on-mismatch.
  3. Checks Rscript is on PATH.
  4. Shells out to coembed.R for the heavy lifting.
  5. Logs the run summary into workflow-log if a run is active.

Coembed is a separate skill (not folded into pseudobulk-construct) so it
can be used standalone — e.g. when a user wants a coembed .rds for
non-Taiji analysis. The output is drop-in for `pseudobulk-construct
--use-existing-clusters`, which closes the loop.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
R_COEMBED = SCRIPT_DIR / "coembed.R"


def _attach_log():
    """Soft-attach to the workflow-log skill's active run. None on any failure."""
    try:
        sibling = (
            Path(__file__).resolve().parent.parent.parent
            / "workflow-log" / "scripts"
        )
        if sibling.exists() and str(sibling) not in sys.path:
            sys.path.insert(0, str(sibling))
        from log import TaijiLog  # type: ignore[import-not-found]
        return TaijiLog.attach_active(work_dir=Path.cwd())
    except Exception:
        return None


def _check_rscript() -> str | None:
    """Return path to Rscript or None if not on PATH."""
    return shutil.which("Rscript")


def _run_rscript(script: Path, *args: str) -> int:
    cmd = ["Rscript", "--no-save", "--no-restore", str(script), *args]
    print(f"[coembed] $ {' '.join(cmd)}", file=sys.stderr)
    proc = subprocess.run(cmd, check=False)
    return proc.returncode


def _macos_memory_budget_gb() -> tuple[float, float] | None:
    """Return (physical_ram_gb, free_swap_gb) on macOS, else None.

    Used both for the pre-flight log ("you have 18 GB RAM + 0.7 GB free
    swap") and for the post-run OOM diagnostic when Rscript exits via
    signal. Cheap — two sysctl calls, no I/O on the .rds files.
    """
    if platform.system() != "Darwin":
        return None
    try:
        ram_b = int(subprocess.check_output(
            ["sysctl", "-n", "hw.memsize"], text=True).strip())
        ram_gb = ram_b / 1024 / 1024 / 1024
    except Exception:
        return None
    swap_free_gb = 0.0
    try:
        line = subprocess.check_output(
            ["sysctl", "-n", "vm.swapusage"], text=True).strip()
        # "total = 24576.00M  used = 23863.81M  free = 712.19M  (encrypted)"
        if "free = " in line:
            v = line.split("free = ", 1)[1].split("M", 1)[0].strip()
            swap_free_gb = float(v) / 1024.0
    except Exception:
        pass
    return (ram_gb, swap_free_gb)


def _diagnose_oom(rc: int) -> str | None:
    """Translate a SIGKILL/SIGTERM/137/143 exit into an OOM mitigation hint.

    The kernel killing R for jetsam pressure does NOT route through R's
    error handler — stderr just stops mid-stage. Without this, the user
    sees `coembed.R exited -15` after 5+ min and has to guess. Return a
    multi-line diagnostic string, or None if the exit looks like a
    normal R error (in which case R's own message is sufficient).
    """
    # subprocess returncode is -signum on POSIX when killed by signal.
    # 128 + signum is the convention if shell wrappers got involved.
    sig = None
    if rc < 0:
        sig = -rc
    elif rc in (137, 143):           # 128 + 9 (SIGKILL) / 128 + 15 (SIGTERM)
        sig = rc - 128
    if sig not in (9, 15):
        return None
    sig_name = "SIGKILL" if sig == 9 else "SIGTERM"
    lines = [
        f"R was killed by {sig_name} mid-pipeline. No R error was emitted —",
        "this is the kernel terminating the process from outside.",
        "",
        "Likely causes (most → least common on Mac):",
        "  1. macOS jetsam OOM (swap budget exhausted during CCA / joint PCA).",
        "  2. SLURM step-mem cgroup limit hit (SBATCH --mem too low).",
        "  3. Manual `kill` from another terminal.",
    ]
    budget = _macos_memory_budget_gb()
    if budget is not None:
        ram_gb, swap_free_gb = budget
        lines.append("")
        lines.append(
            f"This Mac has {ram_gb:.0f} GB physical RAM, "
            f"{swap_free_gb:.1f} GB free swap right now."
        )
        if ram_gb < 24:
            lines.append(
                "  NOTE: 60k+21k cells coembed peaks around ~10-12 GB RSS plus a"
            )
            lines.append(
                "  similar amount of swap pressure. On <24 GB Macs this is at"
            )
            lines.append("  the edge of feasibility.")
    lines.extend([
        "",
        "Mitigations (any one usually unblocks):",
        "  - Run on SLURM: `sbatch --mem=128G --cpus-per-task=8 ...` is reliable",
        "    for >50k merged cells. The skill is unchanged on Linux.",
        "  - Downsample the RNA reference before passing it in:",
        "      rna_ds <- subset(rna, downsample = 20000)",
        "      saveRDS(rna_ds, 'rna_ds.rds')",
        "    Use the downsampled .rds for --rna; ATAC stays full-size since",
        "    only the RNA reference matters for CCA's working set.",
        "  - Lower --n-pcs from 30 to 20 (smaller joint PCA scaled matrix).",
        "  - Free memory: close browsers / quit big GUI apps before running.",
    ])
    return "\n".join(lines)


def _summarize(output_path: Path) -> dict:
    """Read coembed_summary.json next to the output object, if present."""
    summary_path = output_path.parent / "coembed_summary.json"
    if not summary_path.exists():
        return {"error": f"summary file missing at {summary_path}"}
    try:
        return json.loads(summary_path.read_text())
    except Exception as e:
        return {"error": f"could not parse summary: {e}"}


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Co-embed scRNA-seq + scATAC-seq (separate-assay) into a "
            "shared latent space following the Stuart 2019 / Signac "
            "integrate_atac vignette. Output is a single Seurat .rds with "
            "shared PCA/UMAP and de novo clusters, drop-in for "
            "pseudobulk-construct --use-existing-clusters."
        )
    )
    p.add_argument("--rna",  required=True, type=Path,
                   help="Path to scRNA-seq Seurat .rds (RNA assay with raw counts).")
    p.add_argument("--atac", required=True, type=Path,
                   help="Path to scATAC-seq Seurat .rds (ChromatinAssay required).")
    p.add_argument("--output", required=True, type=Path,
                   help="Where to write coembed.rds. Parent dir will hold "
                        "qc/umap.png and coembed_summary.json.")
    p.add_argument("--genome", required=True,
                   help="hg38 / hg19 / mm10 (drives the EnsDb annotation).")
    p.add_argument("--target-cluster-size", type=int, default=200,
                   help="Target mean cluster size for the resolution binary "
                        "search. Default: 200.")
    p.add_argument("--min-cluster-cells", type=int, default=20,
                   help="Drop clusters with fewer than this many cells. "
                        "Default: 20.")
    p.add_argument("--metadata-cols", default=None,
                   help="Comma-separated metadata cols to preserve in the "
                        "merged object.")
    p.add_argument("--n-pcs", type=int, default=30,
                   help="Joint PCA dims for FindNeighbors / RunUMAP. Default: 30.")
    p.add_argument("--lsi-skip-first", type=int, default=1,
                   help="Skip the first N LSI components on ATAC "
                        "(typically depth-correlated). Default: 1.")
    p.add_argument("--resolution", type=float, default=None,
                   help="Force a single resolution; default: scale-aware "
                        "binary search.")
    p.add_argument("--no-cluster", action="store_true",
                   help="Stop after joint UMAP; don't cluster.")
    p.add_argument("--no-plot", action="store_true",
                   help="Skip the QC UMAP rendering.")
    p.add_argument("--reuse-rna-reductions", action="store_true",
                   help="If RNA object already has pca + umap reductions and "
                        "VFs computed, skip Norm/VF/Scale/PCA/UMAP. Saves "
                        "5-10 min on 60k+ cells. Falls back to full "
                        "preprocessing if any prerequisite is missing.")
    p.add_argument("--reuse-atac-reductions", action="store_true",
                   help="If ATAC object already has lsi + umap reductions, "
                        "skip TF-IDF/SVD/UMAP. Saves 3-5 min on 20k+ cells.")
    p.add_argument("--strict-metadata", action="store_true",
                   help="Fail-loud if --metadata-cols values differ between "
                        "RNA and ATAC inputs (e.g. tissue=spleen vs Spleen). "
                        "Default: warn but continue.")
    p.add_argument("--abort-on-memory-risk", action="store_true",
                   help="Refuse to run if the macOS pre-flight estimates "
                        "peak RAM > available RAM + 0.5*free_swap. Default: "
                        "warn but continue. Use on shared Macs where a "
                        "20-min thrash is worse than a fail-fast.")
    p.add_argument("--fragments", type=Path, default=None,
                   help="Path to fragments.tsv.gz. Required when the path "
                        "embedded inside the ATAC ChromatinAssay is stale "
                        "(the typical 'object built on HPC, scp'd to laptop' "
                        "wart). The skill rebuilds the Fragment handle in "
                        "place. Unused if the embedded path resolves on the "
                        "current machine.")
    p.add_argument("--dry-run", action="store_true",
                   help="Validate inputs + print plan; don't run R.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)

    # Existence checks (fast fail).
    for label, path in [("--rna", args.rna), ("--atac", args.atac)]:
        if not path.exists():
            print(f"ERROR: {label} '{path}' does not exist", file=sys.stderr)
            return 2

    # Convention nudge: outputs should live under runs/<name>/, not inside
    # the input data dir. Putting derived artifacts beside immutable inputs
    # makes cleanup hard, breaks per-run reproducibility, and risks future
    # detect-dataset-type runs reclassifying the directory because of the
    # coembed.rds sitting next to the originals.
    try:
        rna_parent  = args.rna.resolve().parent
        atac_parent = args.atac.resolve().parent
        out_parent  = args.output.resolve().parent
        if (out_parent == rna_parent or out_parent == atac_parent
            or rna_parent in out_parent.parents
            or atac_parent in out_parent.parents):
            print(
                f"WARN: --output ({args.output}) is inside an input data "
                f"directory ({rna_parent} / {atac_parent}). Convention is "
                f"to write coembed outputs under runs/<run_name>/coembed/ "
                f"so derived artifacts stay separate from canonical inputs. "
                f"Continuing anyway.",
                file=sys.stderr,
            )
    except (OSError, ValueError):
        pass  # symlink resolution can fail; convention check is advisory

    # Ensure the output path ends with .rds so downstream skills (pseudobulk,
    # detect-dataset-type) can recognise it by extension.
    if args.output.suffix.lower() != ".rds":
        args.output = args.output.with_suffix(args.output.suffix + ".rds")
        print(f"[coembed] NOTE: appending .rds to output path -> {args.output}",
              file=sys.stderr)

    # Output dir setup.
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Dependency check.
    if not _check_rscript():
        if args.dry_run:
            print("NOTE (--dry-run): Rscript not on PATH; "
                  "would fail without --dry-run.", file=sys.stderr)
        else:
            print("ERROR: Rscript not on PATH. On SLURM try "
                  "`module load r/4.3` or activate the conda env.",
                  file=sys.stderr)
            return 3

    if args.dry_run:
        print("[coembed] dry-run plan:", file=sys.stderr)
        print(f"  RNA   : {args.rna}", file=sys.stderr)
        print(f"  ATAC  : {args.atac}", file=sys.stderr)
        print(f"  output: {args.output}", file=sys.stderr)
        print(f"  genome: {args.genome}", file=sys.stderr)
        print(f"  cluster: {'no' if args.no_cluster else 'yes (target ~%d cells/cluster)' % args.target_cluster_size}",
              file=sys.stderr)
        return 0

    # Build R CLI args.
    r_args: list[str] = [
        "--rna", str(args.rna),
        "--atac", str(args.atac),
        "--output", str(args.output),
        "--genome", args.genome,
        "--target-cluster-size", str(args.target_cluster_size),
        "--min-cluster-cells", str(args.min_cluster_cells),
        "--n-pcs", str(args.n_pcs),
        "--lsi-skip-first", str(args.lsi_skip_first),
    ]
    if args.metadata_cols:
        r_args += ["--metadata-cols", args.metadata_cols]
    if args.resolution is not None:
        r_args += ["--resolution", str(args.resolution)]
    if args.no_cluster:
        r_args += ["--no-cluster"]
    if args.no_plot:
        r_args += ["--no-plot"]
    if args.reuse_rna_reductions:
        r_args += ["--reuse-rna-reductions"]
    if args.reuse_atac_reductions:
        r_args += ["--reuse-atac-reductions"]
    if args.strict_metadata:
        r_args += ["--strict-metadata"]
    if args.abort_on_memory_risk:
        r_args += ["--abort-on-memory-risk"]
    if args.fragments is not None:
        if not args.fragments.exists():
            print(f"ERROR: --fragments '{args.fragments}' does not exist",
                  file=sys.stderr)
            return 2
        r_args += ["--fragments", str(args.fragments)]

    # Surface available memory budget BEFORE Rscript runs, so a later OOM
    # diagnostic has a reference point and the user knows up-front if the
    # situation is risky.
    budget = _macos_memory_budget_gb()
    if budget is not None:
        ram_gb, swap_free_gb = budget
        print(f"[coembed] memory budget: {ram_gb:.0f} GB RAM, "
              f"{swap_free_gb:.1f} GB free swap (macOS)", file=sys.stderr)

    rc = _run_rscript(R_COEMBED, *r_args)
    if rc != 0:
        oom_hint = _diagnose_oom(rc)
        if oom_hint is not None:
            print("", file=sys.stderr)
            print("[coembed] " + "-" * 60, file=sys.stderr)
            for line in oom_hint.splitlines():
                print(f"[coembed] {line}" if line else "[coembed]",
                      file=sys.stderr)
            print("[coembed] " + "-" * 60, file=sys.stderr)
        print(f"ERROR: coembed.R exited {rc}", file=sys.stderr)
        return rc

    # Summarize for the user.
    summary = _summarize(args.output)
    if "error" not in summary:
        print(f"\n[coembed] success.", file=sys.stderr)
        print(f"  cells   : {summary.get('n_cells_total')} "
              f"({summary.get('n_cells_rna')} RNA + "
              f"{summary.get('n_cells_atac')} ATAC)", file=sys.stderr)
        print(f"  anchors : {summary.get('n_anchors')}", file=sys.stderr)
        print(f"  clusters: kept {summary.get('n_clusters_kept')} "
              f"of {summary.get('n_clusters_total')} "
              f"(resolution {summary.get('chosen_resolution')})",
              file=sys.stderr)
        print(f"  output  : {args.output}", file=sys.stderr)
        print(f"\n  next step: pseudobulk-construct "
              f"--input {args.output} --use-existing-clusters --fragments <fragments.tsv.gz> "
              f"--genome {args.genome} ...", file=sys.stderr)

    # Best-effort log entry.
    log = _attach_log()
    if log is not None:
        try:
            log.append_custom(
                "coembed_construct",
                f"Stage SC.1 — coembed RNA + ATAC ({args.genome})",
                summary,
                status="pass" if "error" not in summary else "fail",
            )
        except Exception:
            pass

    return 0


if __name__ == "__main__":
    sys.exit(main())
