"""Structured validation of a sample sheet prior to launch.

Returns a `ValidationReport` with `errors` and `warnings`. The CLI exits
non-zero if errors are present. Warnings do not block but are surfaced
prominently so the user can decide whether to proceed with `--force`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from . import detect
from .samplesheet import Sample, SampleSheet


@dataclass
class ValidationReport:
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors

    def extend(self, other: "ValidationReport") -> None:
        self.errors.extend(other.errors)
        self.warnings.extend(other.warnings)


def _check_file(path: Path | None, label: str, sample_id: str, optional: bool = False) -> list[str]:
    if path is None:
        if optional:
            return []
        return [f"sample '{sample_id}': missing {label}"]
    if not path.exists():
        return [f"sample '{sample_id}': {label} does not exist: {path}"]
    if not path.is_file():
        # Directory is OK for sc inputs that point at a 10x-style folder; only
        # complain if it's neither a file nor a directory.
        if not path.is_dir():
            return [f"sample '{sample_id}': {label} is not a file or directory: {path}"]
    return []


def validate_sample(sample: Sample) -> ValidationReport:
    rep = ValidationReport()

    # RNA + ATAC pairing — both sides must be present (this is the core invariant
    # that distinguishes Taiji inputs from single-assay pipelines).
    rep.errors.extend(_check_file(sample.rna_r1, "rna_r1", sample.sample_id))
    rep.errors.extend(_check_file(sample.atac_r1, "atac_r1", sample.sample_id))
    # R2 is optional (single-end bulk is valid) but if present, must exist.
    rep.errors.extend(_check_file(sample.rna_r2, "rna_r2", sample.sample_id, optional=True))
    rep.errors.extend(_check_file(sample.atac_r2, "atac_r2", sample.sample_id, optional=True))

    # Assay-type cross-check (heuristic; warning-only).
    det = detect.check(sample)
    if det.warning:
        rep.warnings.append(det.warning)

    # For single-cell and multiome, single-end reads are a red flag.
    if sample.assay_type in {"sc", "multiome"}:
        if sample.rna_r2 is None and sample.rna_r1.suffix.lower() in {".fastq", ".gz", ".fq"}:
            rep.warnings.append(
                f"sample '{sample.sample_id}': sc/multiome declared but RNA is single-end FASTQ "
                f"— most 10x pipelines require paired reads with cell barcodes in R1."
            )

    return rep


def validate_sheet(sheet: SampleSheet) -> ValidationReport:
    rep = ValidationReport()

    if len(sheet) == 0:
        rep.errors.append("sample sheet has zero samples")
        return rep

    # Genome consistency: Taiji cross-references samples, so mixing builds is
    # almost always an error. We flag as an error, not a warning.
    if len(sheet.genomes) > 1:
        rep.errors.append(
            f"sample sheet mixes genome builds: {sorted(sheet.genomes)}. "
            f"Taiji expects a single reference per run."
        )

    # Assay routing: bulk and sc Taiji have different downstream pipelines.
    # Multiome is compatible with sc. Mixing bulk with sc/multiome in one run
    # is not supported — split into two sheets.
    if {"bulk"} & sheet.assay_types and ({"sc", "multiome"} & sheet.assay_types):
        rep.errors.append(
            f"sample sheet mixes bulk and single-cell assays: {sorted(sheet.assay_types)}. "
            f"Run them as separate Taiji jobs."
        )

    for s in sheet.samples:
        rep.extend(validate_sample(s))

    return rep
