"""Bulk vs single-cell detection.

Policy: the sample sheet's `assay_type` column is authoritative. These heuristics
run as a *cross-check* and only produce warnings. They never silently override
the declared type — doing so would make runs non-reproducible.

Heuristic signals
-----------------
Single-cell ATAC:
  * A `fragments.tsv.gz` file anywhere in the ATAC path's directory
  * A `barcodes.tsv` / `barcodes.tsv.gz`
  * 10x output folder layout (`filtered_feature_bc_matrix/` next to the path)
  * File name contains "fragments" or "atac_peak_bc_matrix"

Single-cell RNA:
  * `.h5ad`, `.loom`, or `.rds` (Seurat/SCE) extension
  * 10x filtered/raw matrix directory
  * File name contains "cellranger", "starsolo", "kb_python"

Multiome: both RNA and ATAC point at a 10x-Multiome-shaped directory, OR file
names contain "multiome".

Bulk: no single-cell signals, inputs look like raw FASTQ.

Limitations
-----------
These heuristics are best-effort. Edge cases this will miss:
  * low-input bulk ATAC that sparsely resembles sc
  * custom pipelines that rename fragment files
  * pre-aligned BAM inputs (we do not currently inspect BAMs)
Treat warnings as prompts to double-check, not as a ground truth.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from .samplesheet import AssayType, Sample

SC_ATAC_FILE_HINTS = ("fragments.tsv", "fragments.tsv.gz", "atac_peak_bc_matrix")
SC_RNA_EXTS = (".h5ad", ".loom", ".rds")
TENX_DIR_HINTS = ("filtered_feature_bc_matrix", "raw_feature_bc_matrix")
MULTIOME_HINTS = ("multiome",)

Verdict = Literal["bulk", "sc", "multiome", "unknown"]


@dataclass
class DetectionResult:
    sample_id: str
    declared: AssayType
    inferred: Verdict
    agrees: bool
    signals: list[str]

    @property
    def warning(self) -> str | None:
        if self.agrees or self.inferred == "unknown":
            return None
        return (
            f"sample '{self.sample_id}': declared assay_type={self.declared!r} but "
            f"heuristics suggest {self.inferred!r} (signals: {', '.join(self.signals)})"
        )


def _path_hints(p: Path) -> list[str]:
    """Return a list of heuristic signals found on or near `p`."""
    hints: list[str] = []
    name = p.name.lower()
    parent = p.parent if p.exists() or p.parent.exists() else Path()

    if p.suffix.lower() in SC_RNA_EXTS:
        hints.append(f"sc-rna-ext:{p.suffix}")
    for hint in SC_ATAC_FILE_HINTS:
        if hint in name:
            hints.append(f"sc-atac-name:{hint}")
    for hint in TENX_DIR_HINTS:
        if hint in name or (parent and hint in parent.name.lower()):
            hints.append(f"tenx-dir:{hint}")
    for hint in MULTIOME_HINTS:
        if hint in name or (parent and hint in parent.name.lower()):
            hints.append(f"multiome-name:{hint}")

    # Sibling files (only if parent dir is readable)
    try:
        if parent.is_dir():
            siblings = {s.name.lower() for s in parent.iterdir()}
            if any("fragments.tsv" in s for s in siblings):
                hints.append("sibling:fragments.tsv")
            if any(s in siblings for s in ("barcodes.tsv", "barcodes.tsv.gz")):
                hints.append("sibling:barcodes.tsv")
    except (PermissionError, OSError):
        pass

    return hints


def infer_assay(sample: Sample) -> tuple[Verdict, list[str]]:
    """Infer an assay type from file layout. Returns (verdict, signals)."""
    rna_hints = _path_hints(sample.rna_r1)
    atac_hints = _path_hints(sample.atac_r1)
    all_hints = rna_hints + atac_hints

    multiome = any("multiome" in h for h in all_hints)
    sc_atac = any(h.startswith("sc-atac-name") or h == "sibling:fragments.tsv" for h in atac_hints)
    sc_rna = any(h.startswith("sc-rna-ext") for h in rna_hints) or any(
        h.startswith("tenx-dir") for h in rna_hints
    )

    if multiome:
        return "multiome", all_hints
    if sc_atac and sc_rna:
        return "sc", all_hints
    if sc_atac or sc_rna:
        # Mixed: one side looks sc, the other doesn't. Flag as unknown so the
        # caller can decide — this is a common human error worth surfacing.
        return "unknown", all_hints
    if not all_hints:
        return "bulk", []
    return "unknown", all_hints


def check(sample: Sample) -> DetectionResult:
    inferred, signals = infer_assay(sample)
    agrees = inferred == sample.assay_type or inferred == "unknown"
    return DetectionResult(
        sample_id=sample.sample_id,
        declared=sample.assay_type,
        inferred=inferred,
        agrees=agrees,
        signals=signals,
    )
