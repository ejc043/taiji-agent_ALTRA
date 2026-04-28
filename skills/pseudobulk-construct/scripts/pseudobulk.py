"""Orchestrator for pseudobulk construction.

Thin Python CLI that:
  1. Gates on detect-dataset-type (refuses bulk / sc-undetermined / separate-assay
     without transferred labels).
  2. Checks that Rscript + MACS2/3 are on PATH.
  3. Shells out to R scripts for object loading, WNN clustering, resolution
     binary search, per-cluster RNA aggregation, and per-cluster ATAC peak
     calling via Signac::CallPeaks.
  4. Emits a manifest.tsv shaped exactly like `build-taiji-input --samples`
     expects, so the two skills chain cleanly.

The heavy lifting lives in sibling R scripts; this file is deliberately short
and does no domain-specific computation itself.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

# R scripts the orchestrator drives (kept as siblings; paths are resolved lazily
# so a user can relocate the skill's scripts/ dir without breaking the CLI).
R_LOAD_AND_CLUSTER = SCRIPT_DIR / "load_and_cluster.R"
R_AGGREGATE_RNA = SCRIPT_DIR / "aggregate_rna.R"
R_CALL_PEAKS = SCRIPT_DIR / "call_peaks.R"


def _attach_log():
    """Soft-attach to the workflow-log skill's active run. Returns None on
    any failure (skill missing, no active run, import error). Never raises."""
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


def _summarize_run(args, gate: dict, groups: list[dict],
                   plan: dict, manifest_path: Path) -> dict:
    """Pull together the dict that gets logged to workflow-log."""
    trace_path = args.output_dir / "resolution_trace.json"
    trace = {}
    if trace_path.exists():
        try:
            trace = json.loads(trace_path.read_text())
        except Exception:
            pass
    sizes_dropped = {}
    # plan.dropped_clusters carries the IDs but not sizes. Sizes are in
    # resolution_trace's last iteration only as a count, so we skip per-id
    # detail unless load_and_cluster.R writes it explicitly later.
    for cl in plan.get("dropped_clusters") or []:
        sizes_dropped[cl] = None  # signal "dropped, size unknown from this view"
    return {
        "clustering_signal": trace.get("signal") or args.clustering_signal,
        "n_cells_in":        trace.get("n_cells"),
        "chosen_resolution": trace.get("chosen_resolution"),
        "resolution_trace":  trace.get("iterations", []),
        "n_clusters_kept":   len(plan.get("kept_clusters") or []),
        "n_clusters_dropped":len(plan.get("dropped_clusters") or []),
        "dropped_cluster_sizes": sizes_dropped,
        "n_groups_planned":  len(groups),
        "manifest_path":     str(manifest_path),
        "metadata_cols":     plan.get("metadata_cols"),
        "transferred_label_col": args.transferred_label_col,
        "coembedding_required": gate.get("sc_modality") == "separate-assay",
        "gate_classification":  gate.get("classification"),
        "gate_sc_modality":     gate.get("sc_modality"),
    }


# ---------------------------------------------------------------------------
# Dependency checks
# ---------------------------------------------------------------------------


@dataclass
class DepCheck:
    ok: bool
    missing: list[str] = field(default_factory=list)
    versions: dict[str, str] = field(default_factory=dict)


def resolve_peak_caller(user_choice: str | None) -> str | None:
    """Pick a peak caller. If user passed an explicit choice, honor it.
    Otherwise auto-detect: prefer macs3 (typical when installed via the
    Signac/Bioconductor stack), fall back to macs2. Returns None if neither
    is on PATH (caller surfaces this as an error)."""
    if user_choice:
        return user_choice if shutil.which(user_choice) else None
    for candidate in ("macs3", "macs2"):
        if shutil.which(candidate):
            return candidate
    return None


def check_dependencies(peak_caller: str | None) -> DepCheck:
    """Check that Rscript and the chosen peak caller are on PATH.

    Does NOT check R package availability — that's the R script's job and it
    emits a clearer error than we could. We only verify the interpreters are
    reachable so the user gets a fast "you forgot to `module load r`" message
    instead of a confusing subprocess failure deep in the run.

    `peak_caller` may be None when auto-detection found neither macs2 nor
    macs3 on PATH; that's reported in `missing` as 'macs2|macs3'.
    """
    check = DepCheck(ok=True)
    binaries: list[str] = ["Rscript"]
    if peak_caller is None:
        check.ok = False
        check.missing.append("macs2|macs3")
    else:
        binaries.append(peak_caller)
    for binary in binaries:
        path = shutil.which(binary)
        if path is None:
            check.ok = False
            check.missing.append(binary)
        else:
            check.versions[binary] = path
    return check


# ---------------------------------------------------------------------------
# detect-dataset-type gate
# ---------------------------------------------------------------------------


def _load_detector():
    """Soft-import the sibling detect-dataset-type skill's library API.

    Returns None if the skill is not installed alongside this one — the user
    may have only installed pseudobulk-construct, in which case we skip the
    gate with a warning rather than hard-failing.
    """
    sibling = (
        Path(__file__).resolve().parent.parent.parent
        / "detect-dataset-type"
        / "scripts"
    )
    if sibling.exists() and str(sibling) not in sys.path:
        sys.path.insert(0, str(sibling))
    try:
        from detect import detect_dataset_type  # type: ignore[import-not-found]

        return detect_dataset_type
    except ImportError:
        return None


def gate_on_detect(input_path: Path, fragments_path: Path | None) -> dict:
    """Run the detect-dataset-type library on the input's parent directory
    (and the fragments dir if supplied) and decide whether to proceed.

    Raises SystemExit(2) with a clear message if the input is not usable.
    Returns a dict with `classification`, `sc_modality`, and `scan_roots` on
    success — the caller uses sc_modality to pick the clustering signal.
    """
    detector = _load_detector()
    if detector is None:
        print(
            "NOTE: detect-dataset-type skill not found; skipping pre-flight. "
            "Pass --skip-data-type-check to silence this notice.",
            file=sys.stderr,
        )
        return {"classification": "unknown", "sc_modality": None, "scan_roots": []}

    roots = [input_path if input_path.is_dir() else input_path.parent]
    if fragments_path is not None:
        frag_root = fragments_path.parent
        if frag_root not in roots:
            roots.append(frag_root)

    result = detector([str(r) for r in roots], recursive=True, strict_mixed=True)

    if result.classification == "bulk":
        print(
            f"ERROR: {input_path} looks like bulk data, not single-cell. "
            "This skill only runs on .rds / .h5ad / .h5mu objects. Use "
            "build-taiji-input directly on the bulk files.",
            file=sys.stderr,
        )
        sys.exit(2)

    if result.classification == "mixed":
        print(
            f"ERROR: {input_path} has both bulk and single-cell signatures. "
            "Split the directories before running pseudobulk-construct.",
            file=sys.stderr,
        )
        sys.exit(2)

    if result.classification == "unknown":
        print(
            f"WARNING: detect-dataset-type could not classify {input_path}. "
            "Proceeding on the assumption that --input is the SC object.",
            file=sys.stderr,
        )

    return {
        "classification": result.classification,
        "sc_modality": getattr(result, "sc_modality", None),
        "scan_roots": [str(r) for r in roots],
    }


def pick_clustering_signal(
    user_choice: str | None,
    sc_modality: str | None,
) -> str:
    """Pick 'wnn' | 'rna' | 'atac' based on sc_modality unless user overrode.

    Refuses (SystemExit 2) on sc-undetermined or separate-assay+unclear-labels
    because those cases need a human decision before we burn R time.
    """
    if user_choice:
        return user_choice

    if sc_modality in ("multiome", None):
        # None means detect-dataset-type wasn't available OR classification
        # was unknown. Default to WNN and let the R script complain if the
        # object is RNA-only.
        return "wnn"

    if sc_modality == "separate-assay":
        # The R script checks for transferred labels explicitly; we pick ATAC
        # LSI as the clustering signal because that's the canonical path
        # post-integrate_atac (cluster on ATAC, use transferred RNA labels for
        # metadata stratification).
        print(
            "NOTE: separate-assay input detected. Clustering on ATAC LSI. "
            "Transferred RNA labels (default column: predicted.id) will be "
            "used as a metadata axis. If you haven't run Signac's "
            "integrate_atac workflow yet, this skill will refuse: "
            "https://stuartlab.org/signac/articles/integrate_atac",
            file=sys.stderr,
        )
        return "atac"

    if sc_modality == "sc-undetermined":
        print(
            "ERROR: detect-dataset-type could not determine whether this is "
            "multiome or separate-assay. Rerun with an explicit "
            "--clustering-signal flag, or disambiguate the input by "
            "renaming/providing a .h5mu. See: "
            "https://stuartlab.org/signac/articles/integrate_atac",
            file=sys.stderr,
        )
        sys.exit(2)

    # Unreachable under current sc_modality enum; keep for forward-compat.
    return "wnn"


# ---------------------------------------------------------------------------
# R + MACS2 invocation helpers
# ---------------------------------------------------------------------------


def run_rscript(script: Path, *args: str) -> None:
    """Run `Rscript SCRIPT ARGS...`, stream stdout/stderr, raise on non-zero."""
    cmd = ["Rscript", "--no-save", "--no-restore", str(script), *args]
    print(f"[pseudobulk] $ {' '.join(cmd)}", file=sys.stderr)
    proc = subprocess.run(cmd, check=False)
    if proc.returncode != 0:
        print(f"ERROR: {script.name} exited with {proc.returncode}", file=sys.stderr)
        sys.exit(proc.returncode)


# (Per-group ATAC peak calling now lives in call_peaks.R via Signac::CallPeaks.
# Python no longer drives MACS2 directly — moving it to R buys us Signac's
# CreateFragmentObject plumbing and per-group barcode reconciliation against
# the fragments file. See call_peaks.R for the suffix-stripping logic.)


# ---------------------------------------------------------------------------
# Manifest writer
# ---------------------------------------------------------------------------


def write_manifest(
    output_dir: Path,
    groups: list[dict],
    genome: str,
) -> Path:
    """Emit a TSV manifest with the exact columns `build-taiji-input --samples`
    expects: name, cohort, group, rna_seq, atac_seq, genome.

    `groups` is a list of dicts produced by the R scripts, one per
    (cluster × metadata_col × metadata_value) triple with the corresponding
    output paths.
    """
    manifest = output_dir / "manifest.tsv"
    fields = ["name", "cohort", "group", "rna_seq", "atac_seq", "genome"]
    with manifest.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for g in groups:
            w.writerow(
                {
                    "name": g["name"],
                    "cohort": g["cohort"],
                    "group": g["group"],
                    "rna_seq": g.get("rna_seq", ""),
                    "atac_seq": g.get("atac_seq", ""),
                    "genome": genome,
                }
            )
    return manifest


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Construct pseudobulk RNA TSVs + ATAC narrowPeak files from a "
            "single-cell object, ready to feed build-taiji-input."
        )
    )
    p.add_argument("--input", required=True, type=Path,
                   help="Path to .rds / .h5ad / .h5mu single-cell object.")
    p.add_argument("--fragments", type=Path, default=None,
                   help="Path to fragments.tsv.gz (+ .tbi). Required for ATAC peaks.")
    p.add_argument("--genome", required=True,
                   help="Genome tag (hg19, hg38, mm9, mm10, mm39).")
    p.add_argument("--output-dir", required=True, type=Path,
                   help="Directory to write outputs.")
    p.add_argument("--metadata-cols", default=None,
                   help="Comma-separated metadata columns to stratify on. "
                        "Multiple columns produce the FULL CROSS-PRODUCT "
                        "within each cluster (e.g. genotype,tissue -> per "
                        "cluster, all of {WT-spleen, WT-siiel, KO-spleen, "
                        "KO-siiel}). Default: auto-detect.")
    p.add_argument("--cohort-col", default=None,
                   help="Which --metadata-cols column's value becomes the "
                        "manifest 'cohort' label (the axis Taiji compares "
                        "across). Default: first --metadata-cols entry.")
    p.add_argument("--target-cluster-size", type=int, default=200,
                   help="Target mean cluster size for resolution search. "
                        "Default: 200.")
    p.add_argument("--min-cluster-cells", type=int, default=20,
                   help="Drop clusters with fewer than this many cells. Default: 20.")
    p.add_argument("--clustering-signal", choices=["wnn", "rna", "atac"],
                   default=None,
                   help="Override auto-selected signal (based on sc_modality).")
    p.add_argument("--transferred-label-col", default="predicted.id",
                   help="meta.data column holding Signac-transferred RNA labels. "
                        "Default: predicted.id.")
    p.add_argument("--peak-caller", choices=["macs2", "macs3"], default=None,
                   help="Peak caller binary. Default: auto-detect (prefers "
                        "macs3, then macs2, whichever is on PATH).")
    p.add_argument("--rna-only", action="store_true",
                   help="Skip ATAC peak calling even if fragments are supplied.")
    p.add_argument("--atac-only", action="store_true",
                   help="Skip RNA aggregation.")
    p.add_argument("--skip-data-type-check", action="store_true",
                   help="Bypass the detect-dataset-type pre-flight.")
    p.add_argument("--yes", action="store_true",
                   help="Non-interactive: accept auto-detected metadata columns "
                        "without prompting.")
    p.add_argument("--dry-run", action="store_true",
                   help="Plan only: validate inputs and print the resolved "
                        "plan (no R/MACS invocations). Useful for checking "
                        "the wiring before submitting a long SLURM job.")
    p.add_argument("--no-plot", action="store_true",
                   help="Skip the QC UMAP rendering in load_and_cluster.R "
                        "(default: write qc/umap.png with panels for clusters, "
                        "assay [if present], and detected metadata cols).")
    p.add_argument("--use-existing-clusters", action="store_true",
                   help="Skip clustering — use the seurat_clusters column "
                        "already in the input object (e.g. when the input is "
                        "a coembed-construct output). Without this flag, the "
                        "skill re-clusters on the chosen --clustering-signal, "
                        "which discards any pre-computed shared-space "
                        "clustering you already did.")
    args = p.parse_args(argv)

    if args.rna_only and args.atac_only:
        p.error("--rna-only and --atac-only are mutually exclusive.")
    return args


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)

    # 1. detect-dataset-type gate. Two reasons this is conditional now:
    #    - --skip-data-type-check: explicit user opt-out (preserved).
    #    - --use-existing-clusters: the user is asserting the input is a
    #      pre-clustered SC object (e.g. coembed-construct's output). The
    #      gate's modality auto-routing is moot in that case, and the gate
    #      can produce confusing reclassifications when the coembed.rds
    #      lives next to the original RNA + ATAC .rds files (it'd see
    #      "separate-assay" again and try to re-pick a clustering signal).
    #      We still want to refuse hard on bulk inputs, so do a single-file
    #      classification rather than scanning the whole parent directory.
    if args.skip_data_type_check:
        gate = {"classification": "skipped", "sc_modality": None}
    elif args.use_existing_clusters:
        # Lightweight gate: just confirm the input file extension is one
        # of the SC formats. Don't scan the parent dir.
        ext = args.input.suffix.lower()
        if ext in (".rds", ".h5ad", ".h5mu"):
            gate = {"classification": "single-cell",
                    "sc_modality": "pre-coembedded"}
        else:
            print(f"ERROR: --use-existing-clusters requires an SC object "
                  f"(.rds / .h5ad / .h5mu); got '{ext}'.", file=sys.stderr)
            return 2
    else:
        gate = gate_on_detect(args.input, args.fragments)

    # Pick the clustering signal. With --use-existing-clusters the R script
    # ignores the signal value (it short-circuits to "use the existing
    # seurat_clusters"), so we just pick a benign default to satisfy the
    # CLI contract.
    if args.use_existing_clusters:
        signal = args.clustering_signal or "rna"
    else:
        signal = pick_clustering_signal(args.clustering_signal,
                                        gate.get("sc_modality"))

    # 2. Resolve peak caller (auto-detect macs3 then macs2 unless user gave
    #    an explicit --peak-caller). This is cheap and useful to print even
    #    in --dry-run.
    resolved_peak_caller = resolve_peak_caller(args.peak_caller)
    if resolved_peak_caller and not args.peak_caller:
        print(f"[pseudobulk] auto-selected peak caller: {resolved_peak_caller}",
              file=sys.stderr)
    args.peak_caller = resolved_peak_caller or args.peak_caller or "macs3"

    # 3. --dry-run short-circuit. Print the resolved plan and exit cleanly
    #    WITHOUT invoking any R script. Real planning, not "we tried to run
    #    and it crashed."
    if args.dry_run:
        print("\n=== pseudobulk-construct --dry-run plan ===", file=sys.stderr)
        print(f"  input               : {args.input}", file=sys.stderr)
        print(f"  fragments           : {args.fragments}", file=sys.stderr)
        print(f"  output_dir          : {args.output_dir}", file=sys.stderr)
        print(f"  genome              : {args.genome}", file=sys.stderr)
        print(f"  gate classification : {gate.get('classification')}",
              file=sys.stderr)
        print(f"  gate sc_modality    : {gate.get('sc_modality')}",
              file=sys.stderr)
        print(f"  clustering signal   : {signal}", file=sys.stderr)
        print(f"  use-existing-clusters: {args.use_existing_clusters}",
              file=sys.stderr)
        print(f"  metadata cols       : {args.metadata_cols or '(auto-detect)'}",
              file=sys.stderr)
        print(f"  cohort col          : {args.cohort_col or '(first metadata col)'}",
              file=sys.stderr)
        print(f"  peak caller         : {args.peak_caller}", file=sys.stderr)
        print(f"  rna-only / atac-only: {args.rna_only} / {args.atac_only}",
              file=sys.stderr)
        print("[pseudobulk] --dry-run: skipping all R + MACS invocations.",
              file=sys.stderr)
        return 0

    # 4. Dependency check (real run only — dry-run already returned).
    dep = check_dependencies(resolved_peak_caller)
    if not dep.ok:
        print(
            "ERROR: required binaries missing from PATH: "
            f"{dep.missing}. On SLURM, try `module load r/4.3 "
            "macs2/2.2.9.1` or activate the appropriate conda env. "
            "Use --dry-run to plan without running R/MACS2.",
            file=sys.stderr,
        )
        return 3

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # 3. Load the object, WNN/RNA/ATAC cluster, scale-aware resolution search,
    #    write clusters.csv + resolution_trace.json + per_cluster_barcodes/.
    run_rscript(
        R_LOAD_AND_CLUSTER,
        "--input", str(args.input),
        "--output-dir", str(args.output_dir),
        "--signal", signal,
        "--target-cluster-size", str(args.target_cluster_size),
        "--min-cluster-cells", str(args.min_cluster_cells),
        "--transferred-label-col", args.transferred_label_col,
        *( ["--metadata-cols", args.metadata_cols] if args.metadata_cols else [] ),
        *( ["--cohort-col", args.cohort_col] if args.cohort_col else [] ),
        *( ["--use-existing-clusters"] if args.use_existing_clusters else [] ),
        *( ["--yes"] if args.yes else [] ),
        *( ["--no-plot"] if args.no_plot else [] ),
    )

    # 4. Read the per-group plan the R script wrote.
    plan_path = args.output_dir / "groups_plan.json"
    if not plan_path.exists():
        print(
            f"ERROR: {plan_path} missing; the R clustering step did not "
            "produce a plan. Check Rscript stderr above.",
            file=sys.stderr,
        )
        return 4
    with plan_path.open() as fh:
        plan = json.load(fh)

    groups: list[dict] = plan["groups"]

    # 5. RNA aggregation (R-side): writes rna/<name>.tsv per group.
    if not args.atac_only:
        run_rscript(
            R_AGGREGATE_RNA,
            "--input", str(args.input),
            "--clusters", str(args.output_dir / "clusters.csv"),
            "--groups", str(plan_path),
            "--output-dir", str(args.output_dir / "rna"),
        )
        for g in groups:
            g["rna_seq"] = str(args.output_dir / "rna" / f"{g['name']}.tsv")

    # 6. ATAC peak calling (R-side, via Signac::CallPeaks). One Rscript call
    #    handles all groups; per-group barcode reconciliation against the
    #    fragments file (suffix-stripping), Fragment-object plumbing, and
    #    narrowPeak emission all happen there.
    if not args.rna_only and args.fragments is not None:
        atac_dir = args.output_dir / "atac"
        atac_dir.mkdir(exist_ok=True)
        macs_path = shutil.which(args.peak_caller)
        if macs_path is None:
            print(
                f"ERROR: --peak-caller '{args.peak_caller}' not on PATH; "
                "cannot run Signac::CallPeaks.",
                file=sys.stderr,
            )
            return 5
        run_rscript(
            R_CALL_PEAKS,
            "--input",       str(args.input),
            "--clusters",    str(args.output_dir / "clusters.csv"),
            "--groups",      str(plan_path),
            "--fragments",   str(args.fragments),
            "--output-dir",  str(atac_dir),
            "--genome",      args.genome,
            "--macs-path",   macs_path,
        )
        for g in groups:
            narrow = atac_dir / f"{g['name']}.narrowPeak"
            if narrow.exists() and narrow.stat().st_size > 0:
                g["atac_seq"] = str(narrow)

    # 7. Write the build-taiji-input manifest.
    manifest = write_manifest(args.output_dir, groups, args.genome)
    print(f"[pseudobulk] manifest written: {manifest}", file=sys.stderr)
    print(
        "Next step: run build-taiji-input with "
        f"--samples {manifest} --data-dir {args.output_dir} "
        f"--genome {args.genome}",
        file=sys.stderr,
    )

    # 8. Best-effort log entry.
    log = _attach_log()
    if log is not None:
        try:
            log.append_pseudobulk(_summarize_run(args, gate, groups, plan, manifest))
        except Exception:
            pass

    return 0


if __name__ == "__main__":
    sys.exit(main())
