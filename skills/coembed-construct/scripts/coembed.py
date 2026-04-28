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
                   help="hg38 / hg19 / mm10 / mm39 (drives the EnsDb annotation).")
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

    rc = _run_rscript(R_COEMBED, *r_args)
    if rc != 0:
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
