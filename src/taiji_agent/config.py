"""Generate a Taiji YAML config from a validated sample sheet.

Taiji's YAML format differs between bulk and single-cell flavors and has
evolved between versions. This module emits a *skeleton* that covers the
fields most users customize; site-specific fields (genome index paths, motif
DB, output_dir) are injected from a site config passed on the command line.

IMPORTANT: the exact schema Taiji expects depends on the installed Taiji
version. Treat the emitted YAML as a starting point and validate it with
`taiji run --dry-run` before committing it to a production pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from .samplesheet import SampleSheet


@dataclass
class SiteConfig:
    """Cluster- and lab-specific settings injected into the Taiji config."""

    output_dir: Path
    genome_index: dict[str, Path]      # {"hg38": Path(...), "mm10": Path(...)}
    motif_db: Path
    taiji_binary: str = "taiji"        # or a module-loaded path on HPC
    threads_per_sample: int = 8

    @classmethod
    def from_yaml(cls, path: str | Path) -> "SiteConfig":
        with open(path) as fh:
            data = yaml.safe_load(fh)
        return cls(
            output_dir=Path(data["output_dir"]),
            genome_index={k: Path(v) for k, v in data["genome_index"].items()},
            motif_db=Path(data["motif_db"]),
            taiji_binary=data.get("taiji_binary", "taiji"),
            threads_per_sample=int(data.get("threads_per_sample", 8)),
        )


def _sample_entry(sample, site: SiteConfig) -> dict[str, Any]:
    entry: dict[str, Any] = {
        "id": sample.sample_id,
        "group": sample.condition,
        "rna-seq": {"fastq": [str(sample.rna_r1)]},
        "atac-seq": {"fastq": [str(sample.atac_r1)]},
    }
    if sample.rna_r2:
        entry["rna-seq"]["fastq"].append(str(sample.rna_r2))
    if sample.atac_r2:
        entry["atac-seq"]["fastq"].append(str(sample.atac_r2))
    entry.update(sample.extras)
    return entry


def generate_config(sheet: SampleSheet, site: SiteConfig) -> dict[str, Any]:
    if len(sheet.genomes) != 1:
        raise ValueError("generate_config expects a single-genome sample sheet")
    (genome,) = sheet.genomes

    if genome not in site.genome_index:
        raise ValueError(
            f"site config has no genome_index entry for {genome!r}. "
            f"Available: {sorted(site.genome_index)}"
        )

    # This shape is a common-denominator skeleton; adapt to the specific Taiji
    # version in use. See docstring.
    return {
        "output_dir": str(site.output_dir),
        "genome": genome,
        "genome_index": str(site.genome_index[genome]),
        "motif_file": str(site.motif_db),
        "threads": site.threads_per_sample,
        "input": [_sample_entry(s, site) for s in sheet.samples],
    }


def write_config(sheet: SampleSheet, site: SiteConfig, out: Path) -> Path:
    cfg = generate_config(sheet, site)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)
    return out
