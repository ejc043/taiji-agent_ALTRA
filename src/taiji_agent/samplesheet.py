"""Sample sheet schema + loader.

The sample sheet is the single source of truth for a Taiji run. It is a CSV/TSV
file with one row per biological sample. RNA-seq and ATAC-seq reads live in the
same row (they describe the same sample, measured two ways).

Required columns
----------------
sample_id       Unique per sheet; stable across reruns.
condition       Free-text group label (used for Taiji contrasts).
genome          Genome build identifier (e.g. hg38, mm10). Must be consistent
                within a run.
assay_type      One of {"bulk", "sc", "multiome"}. See `detect.py` for what
                "multiome" implies (10x Multiome: RNA + ATAC from same cells).
rna_r1          Path to RNA FASTQ R1 (or per-cell matrix for sc/multiome).
rna_r2          Path to RNA FASTQ R2. Optional for single-end bulk.
atac_r1         Path to ATAC FASTQ R1 (or fragments.tsv.gz for sc/multiome).
atac_r2         Path to ATAC FASTQ R2. Optional for single-end bulk.

Optional columns are preserved and passed through to the Taiji config.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterator, Literal

from pydantic import BaseModel, Field, field_validator

AssayType = Literal["bulk", "sc", "multiome"]
REQUIRED_COLUMNS = {
    "sample_id",
    "condition",
    "genome",
    "assay_type",
    "rna_r1",
    "atac_r1",
}


class Sample(BaseModel):
    """One biological sample with paired RNA-seq and ATAC-seq inputs."""

    sample_id: str = Field(..., min_length=1)
    condition: str = Field(..., min_length=1)
    genome: str = Field(..., min_length=1)
    assay_type: AssayType
    rna_r1: Path
    rna_r2: Path | None = None
    atac_r1: Path
    atac_r2: Path | None = None
    # Freeform extras that the Taiji config may reference.
    extras: dict[str, str] = Field(default_factory=dict)

    @field_validator("sample_id")
    @classmethod
    def _no_whitespace(cls, v: str) -> str:
        if any(c.isspace() for c in v):
            raise ValueError("sample_id must not contain whitespace")
        return v


class SampleSheet(BaseModel):
    """Parsed sample sheet. Construct via `SampleSheet.load(path)`."""

    samples: list[Sample]
    source: Path | None = None

    def __iter__(self) -> Iterator[Sample]:  # type: ignore[override]
        return iter(self.samples)

    def __len__(self) -> int:
        return len(self.samples)

    def by_assay(self, assay: AssayType) -> list[Sample]:
        return [s for s in self.samples if s.assay_type == assay]

    @property
    def genomes(self) -> set[str]:
        return {s.genome for s in self.samples}

    @property
    def assay_types(self) -> set[AssayType]:
        return {s.assay_type for s in self.samples}

    @classmethod
    def load(cls, path: str | Path) -> "SampleSheet":
        p = Path(path)
        delim = "\t" if p.suffix.lower() in {".tsv", ".txt"} else ","
        with p.open(newline="") as fh:
            reader = csv.DictReader(fh, delimiter=delim)
            fieldnames = reader.fieldnames or []
            missing = REQUIRED_COLUMNS - set(fieldnames)
            if missing:
                raise ValueError(
                    f"Sample sheet {p} is missing required columns: "
                    f"{sorted(missing)}. Required: {sorted(REQUIRED_COLUMNS)}"
                )
            samples: list[Sample] = []
            for row_idx, row in enumerate(reader, start=2):  # header is row 1
                extras = {
                    k: v for k, v in row.items()
                    if k not in REQUIRED_COLUMNS
                    and k not in {"rna_r2", "atac_r2"}
                    and v not in (None, "")
                }
                try:
                    samples.append(
                        Sample(
                            sample_id=row["sample_id"],
                            condition=row["condition"],
                            genome=row["genome"],
                            assay_type=row["assay_type"],  # type: ignore[arg-type]
                            rna_r1=Path(row["rna_r1"]),
                            rna_r2=Path(row["rna_r2"]) if row.get("rna_r2") else None,
                            atac_r1=Path(row["atac_r1"]),
                            atac_r2=Path(row["atac_r2"]) if row.get("atac_r2") else None,
                            extras=extras,
                        )
                    )
                except Exception as e:
                    raise ValueError(f"{p}: row {row_idx}: {e}") from e

        ids = [s.sample_id for s in samples]
        dupes = {sid for sid in ids if ids.count(sid) > 1}
        if dupes:
            raise ValueError(f"Duplicate sample_id in {p}: {sorted(dupes)}")

        return cls(samples=samples, source=p)
