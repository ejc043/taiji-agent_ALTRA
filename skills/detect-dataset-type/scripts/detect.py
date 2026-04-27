"""Classify a dataset directory as bulk / single-cell / mixed / unknown.

Decision is made purely on file extensions (fast, deterministic, no file opens).
Bulk signatures:       .tsv (RNA-seq GeneQuant), .narrowPeak (ATAC-seq),
                       .bedpe (HiC chromatin loops).
Single-cell signatures: .h5ad (AnnData), .rds (R serialized, usually Seurat/SCE).

Exposes:
- detect_dataset_type(paths, recursive=True, strict_mixed=True) -> DetectResult
  (library API for build-taiji-input and other callers)
- CLI: `python -m detect_dataset_type PATH [PATH ...]` with --format, --no-recursive,
  --allow-mixed.

Exit codes:
  0 = classified (bulk, single-cell, or unknown)
  2 = mixed under strict_mixed=True
  3 = path error (missing / unreadable)
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable

BULK_EXTS = {".tsv", ".narrowpeak", ".bedpe"}
SC_EXTS = {".h5ad", ".rds", ".h5mu"}  # .h5mu = MuData (multi-modal, same cells)

# Filenames that signal single-cell regardless of their extension. These take
# priority over BULK_EXTS so that 10x cellranger / cellranger-arc outputs
# (which contain atac_fragments.tsv.gz, barcodes.tsv, features.tsv, etc.) do
# not get misclassified as bulk RNA-seq because of the `.tsv` suffix.
# Match on case-insensitive substring of basename.
SC_FILENAME_SUBSTRINGS = (
    "fragments.tsv",                   # scATAC fragments (10x/Signac)
    "barcodes.tsv",                    # 10x GEX/ATAC barcodes
    "features.tsv",                    # 10x GEX features
    "genes.tsv",                       # older 10x GEX genes file
    "matrix.mtx",                      # sparse count matrix (10x, STARsolo)
    "filtered_feature_bc_matrix.h5",   # 10x GEX or Multiome combined output
    "filtered_peak_bc_matrix.h5",      # 10x ATAC peak-barcode matrix
    "raw_feature_bc_matrix.h5",        # 10x raw counts
)

# Filename patterns that identify which modality a file covers. Used for
# sub-classification of single-cell data into multiome (RNA + ATAC from the
# same cells) vs separate-assay (need co-embedding via the Signac workflow).
RNA_HINTS = ("rna", "gex", "expression", "scrna", "gene_expression")
ATAC_HINTS = ("atac", "scatac", "peak", "chromatin_accessibility")
MULTIOME_HINTS = ("multiome", "multi_omic", "multi-omic", "multimodal", "arc")

IGNORE_DIRS = {".git", "__pycache__", ".ipynb_checkpoints", "node_modules", ".venv", ".tox"}
IGNORE_NAMES = {".DS_Store", "Thumbs.db"}


@dataclass
class DetectResult:
    classification: str  # "bulk" | "single-cell" | "mixed" | "unknown"
    paths_scanned: list[str] = field(default_factory=list)
    bulk_files: dict[str, list[str]] = field(default_factory=dict)
    sc_files: dict[str, list[str]] = field(default_factory=dict)
    other_exts: dict[str, int] = field(default_factory=dict)
    sizes: dict[str, int] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    # Sub-classification (set only when classification == "single-cell"):
    #   "multiome"        -> RNA + ATAC measured from the SAME cells
    #   "separate-assay"  -> RNA and ATAC from DIFFERENT cell populations; co-embedding needed
    #   "sc-undetermined" -> single-cell data but modality layout cannot be inferred from filenames
    sc_modality: str | None = None

    @property
    def ok_for_bulk_taiji(self) -> bool:
        """Gate used by build-taiji-input. True only for clean bulk datasets."""
        return self.classification == "bulk" and not self.errors

    def summary(self) -> str:
        lines = [f"Classification: {self.classification}"]
        if self.sc_modality:
            lines.append(f"SC modality: {self.sc_modality}")
        lines.append(f"Scanned: {', '.join(self.paths_scanned) or '(none)'}")
        total_bulk = sum(len(v) for v in self.bulk_files.values())
        total_sc = sum(len(v) for v in self.sc_files.values())
        if total_bulk:
            lines.append(f"Bulk files: {total_bulk}")
            for ext, files in sorted(self.bulk_files.items()):
                lines.append(f"  {ext}: {len(files)}")
        if total_sc:
            lines.append(f"Single-cell files: {total_sc}")
            for ext, files in sorted(self.sc_files.items()):
                lines.append(f"  {ext}: {len(files)}")
        if self.errors:
            lines.append("Errors:")
            lines.extend(f"  - {e}" for e in self.errors)
        if self.warnings:
            lines.append("Warnings:")
            lines.extend(f"  - {w}" for w in self.warnings)
        return "\n".join(lines)


def _normalized_ext(path: Path) -> str:
    """Extension, lowercased, with .gz peeled off so .tsv.gz -> .tsv."""
    suffixes = path.suffixes
    if suffixes and suffixes[-1].lower() == ".gz" and len(suffixes) >= 2:
        return suffixes[-2].lower()
    return path.suffix.lower()


def _matches_sc_filename(name: str) -> str | None:
    """Return the matched SC filename substring, or None if nothing matches.

    Matching is case-insensitive against the basename. The substring itself is
    used as the 'extension' key in sc_files so callers can tell, e.g., that
    the match was via `fragments.tsv` rather than via `.h5ad`.
    """
    lower = name.lower()
    for sub in SC_FILENAME_SUBSTRINGS:
        if sub in lower:
            return sub
    return None


def _check_bulk_completeness(result: "DetectResult") -> None:
    """Warn for each missing required bulk modality (RNA, ATAC) and note absent HiC."""
    has_rna = ".tsv" in result.bulk_files
    has_atac = ".narrowpeak" in result.bulk_files
    has_hic = ".bedpe" in result.bulk_files

    if not has_rna:
        result.warnings.append(
            "MISSING REQUIRED: no RNA-seq GeneQuant files (.tsv) found. "
            "Taiji requires one 2-column TSV per sample (no header): "
            "gene_symbol<TAB>expression_value. "
            "These are typically produced by RSEM, STARsolo, or featureCounts."
        )
    if not has_atac:
        result.warnings.append(
            "MISSING REQUIRED: no ATAC-seq peak files (.narrowPeak) found. "
            "Taiji requires ENCODE BED6+4 narrowPeak format — standard MACS2/MACS3 output "
            "from pre-aligned BAMs. Columns: chr, start, end, name, score, strand, "
            "signalValue, pValue, qValue, summit."
        )
    if not has_hic:
        result.warnings.append(
            "no HiC chromatin-loop files (.bedpe) found. HiC is optional but improves "
            "TF→target edge accuracy. Format: 6-column BEDPE "
            "(chr1 start1 end1 chr2 start2 end2, tab-separated, no header)."
        )


def _check_sc_completeness(result: "DetectResult", all_basenames: list[str]) -> None:
    """Warn if no SC object file or no fragments file found for an SC dataset."""
    sc_object_exts = {".rds", ".h5ad", ".h5mu"}
    has_object = any(ext in result.sc_files for ext in sc_object_exts)
    if not has_object:
        result.warnings.append(
            "MISSING REQUIRED: no single-cell object file (.rds / .h5ad / .h5mu) found. "
            "pseudobulk-construct requires a Seurat (.rds), AnnData (.h5ad), "
            "or MuData (.h5mu) object to load cells and run clustering."
        )

    has_fragments = "fragments.tsv" in result.sc_files
    if not has_fragments:
        result.warnings.append(
            "MISSING RECOMMENDED: no ATAC fragments file (fragments.tsv.gz) found. "
            "pseudobulk-construct needs a bgzipped, Tabix-indexed 10x-compatible "
            "fragments file for per-cluster ATAC peak calling. Without it, "
            "only RNA pseudobulks can be generated (use --rna-only)."
        )


def _classify_sc_modality(
    sc_files: dict[str, list[str]],
    all_basenames: list[str],
) -> str:
    """Decide whether SC data is multiome, separate-assay, or undetermined.

    Tier logic (first match wins):
      1. Any `.h5mu` file present                       -> multiome
      2. cellranger-arc signature: filtered_feature_bc_matrix.h5
         AND atac_fragments.tsv both present            -> multiome
      3. A filename contains any MULTIOME_HINTS token   -> multiome
      4. Filenames collectively contain BOTH an RNA_HINTS
         token AND an ATAC_HINTS token                  -> separate-assay
      5. Otherwise                                      -> sc-undetermined

    Notes / limitations:
      - This is a filename-only heuristic; it cannot see inside .h5ad/.rds to
        check for a paired ATAC assay slot. A .h5ad named
        `pbmc_multiome.h5ad` lands in tier 3; a generic `pbmc.h5ad` with no
        ATAC siblings will end up `sc-undetermined`.
      - Tier 4 triggers on any filename with both hints — e.g. `rna.h5ad` +
        `atac.h5ad` in the same tree, which is the common "separate-assay"
        layout that needs Signac/Seurat co-embedding.
    """
    lowered = [n.lower() for n in all_basenames]
    joined = " ".join(lowered)

    # Tier 1: MuData
    if ".h5mu" in sc_files:
        return "multiome"

    # Tier 2: cellranger-arc directory signature
    has_feature_h5 = any("filtered_feature_bc_matrix.h5" in n for n in lowered)
    has_atac_fragments = any("atac_fragments.tsv" in n for n in lowered)
    if has_feature_h5 and has_atac_fragments:
        return "multiome"

    # Tier 3: explicit multiome tokens in filenames
    if any(hint in joined for hint in MULTIOME_HINTS):
        return "multiome"

    # Tier 4: both RNA and ATAC hints -> separate-assay (co-embedding required)
    has_rna = any(hint in joined for hint in RNA_HINTS)
    has_atac = any(hint in joined for hint in ATAC_HINTS)
    if has_rna and has_atac:
        return "separate-assay"

    return "sc-undetermined"


def _iter_files(root: Path, recursive: bool) -> Iterable[Path]:
    if root.is_file():
        yield root
        return
    if not root.is_dir():
        return
    walker = root.rglob("*") if recursive else root.iterdir()
    for p in walker:
        parts = p.parts
        if any(part in IGNORE_DIRS for part in parts):
            continue
        if p.name in IGNORE_NAMES:
            continue
        if p.name.startswith("."):
            continue
        if p.is_file():
            yield p


def detect_dataset_type(
    paths: Iterable[str | Path],
    recursive: bool = True,
    strict_mixed: bool = True,
) -> DetectResult:
    """Classify a set of paths as bulk / single-cell / mixed / unknown.

    Args:
        paths: One or more files or directories to scan.
        recursive: Recurse into subdirectories. Default True.
        strict_mixed: If both bulk and SC signatures appear, mark as an error
            (classification='mixed'). Set False to downgrade to a warning.

    Returns:
        DetectResult. Never raises for classification reasons; path errors go
        to result.errors and `paths_scanned` still reflects what was attempted.
    """
    result = DetectResult(classification="unknown")
    roots = [Path(p).expanduser() for p in paths]
    all_basenames: list[str] = []

    for root in roots:
        result.paths_scanned.append(str(root))
        if not root.exists():
            result.errors.append(f"path does not exist: {root}")
            continue

        for f in _iter_files(root, recursive):
            ext = _normalized_ext(f)
            all_basenames.append(f.name)
            try:
                result.sizes[str(f)] = f.stat().st_size
            except OSError:
                pass

            # SC filename patterns take priority over BULK_EXTS. This prevents
            # 10x cellranger-arc outputs (atac_fragments.tsv.gz, barcodes.tsv,
            # features.tsv, matrix.mtx, *.h5) from being miscounted as bulk.
            sc_match = _matches_sc_filename(f.name)
            if sc_match is not None:
                result.sc_files.setdefault(sc_match, []).append(str(f))
            elif ext in BULK_EXTS:
                result.bulk_files.setdefault(ext, []).append(str(f))
            elif ext in SC_EXTS:
                result.sc_files.setdefault(ext, []).append(str(f))
            else:
                result.other_exts[ext] = result.other_exts.get(ext, 0) + 1

    has_bulk = bool(result.bulk_files)
    has_sc = bool(result.sc_files)

    if has_bulk and has_sc:
        result.classification = "mixed"
        bulk_exts = sorted(result.bulk_files.keys())
        sc_exts = sorted(result.sc_files.keys())
        msg = (
            f"both bulk ({bulk_exts}) and single-cell ({sc_exts}) "
            "signatures present. Bulk Taiji cannot consume mixed datasets — "
            "split the directories or remove the minority files."
        )
        (result.errors if strict_mixed else result.warnings).append(msg)
    elif has_sc:
        result.classification = "single-cell"
        result.sc_modality = _classify_sc_modality(result.sc_files, all_basenames)
        _check_sc_completeness(result, all_basenames)

        if result.sc_modality == "multiome":
            result.warnings.append(
                "single-cell multiome detected (RNA + ATAC from the SAME cells). "
                "Not compatible with build-taiji-input; route to a single-cell "
                "Taiji workflow that accepts paired RNA+ATAC (e.g. sc-Taiji with "
                "a MuData/AnnData + fragments input)."
            )
        elif result.sc_modality == "separate-assay":
            result.warnings.append(
                "single-cell RNA-seq and scATAC-seq detected from DIFFERENT cell "
                "populations. These assays must be co-embedded before downstream "
                "regulatory-network analysis. Use the Signac/Seurat integration "
                "workflow to transfer RNA labels/expression onto the ATAC cells "
                "via reciprocal-LSI anchors: "
                "https://stuartlab.org/signac/articles/integrate_atac . "
                "After co-embedding, export a paired object and feed it to the "
                "single-cell Taiji workflow (build-taiji-input is bulk-only)."
            )
        else:  # sc-undetermined
            result.warnings.append(
                "single-cell data detected but the RNA/ATAC modality layout "
                "could not be inferred from filenames. If this is multiome "
                "(same cells), rename files with a 'multiome' token or provide "
                "a .h5mu. If RNA and ATAC come from DIFFERENT cells, co-embed "
                "via Signac's integrate_atac workflow before downstream "
                "analysis: https://stuartlab.org/signac/articles/integrate_atac . "
                "build-taiji-input is bulk-only and will not run on this dataset."
            )
    elif has_bulk:
        result.classification = "bulk"
        _check_bulk_completeness(result)
        # .tsv alone is a weak signal — could be any tab-separated file.
        if set(result.bulk_files.keys()) == {".tsv"}:
            result.warnings.append(
                ".tsv is the only bulk signature present and is ambiguous. "
                "Confirm these are RNA-seq GeneQuant files (gene IDs in column 1)."
            )
    else:
        result.classification = "unknown"
        top = sorted(result.other_exts.items(), key=lambda kv: -kv[1])[:5]
        result.warnings.append(
            "no recognized bulk or single-cell signatures found. "
            f"Top extensions seen: {top or '(none)'}."
        )

    return result


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Classify a dataset as bulk / single-cell / mixed / unknown "
            "based on file extensions."
        )
    )
    p.add_argument("paths", nargs="+", help="Directory or file path(s) to classify.")
    p.add_argument(
        "--no-recursive",
        action="store_true",
        help="Do not recurse into subdirectories.",
    )
    p.add_argument(
        "--allow-mixed",
        action="store_true",
        help="Report mixed datasets as a warning instead of erroring.",
    )
    p.add_argument(
        "--format",
        choices=["text", "json"],
        default="text",
        help="Output format. Default: text.",
    )
    return p.parse_args(argv)


def _attach_log():
    """Soft-import workflow-log skill. Returns a TaijiLog bound to the active
    run, or None if (a) the skill isn't installed alongside or (b) no active
    run is registered. Never raises — logging never breaks the workflow."""
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


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    result = detect_dataset_type(
        args.paths,
        recursive=not args.no_recursive,
        strict_mixed=not args.allow_mixed,
    )

    if args.format == "json":
        print(json.dumps(asdict(result), indent=2))
    else:
        print(result.summary())

    # Best-effort log entry. Logging failures are swallowed.
    log = _attach_log()
    if log is not None:
        try:
            log.append_detect(result)
        except Exception:
            pass

    # path errors take precedence over mixed
    if any("does not exist" in e for e in result.errors):
        return 3
    if result.classification == "mixed" and result.errors:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
