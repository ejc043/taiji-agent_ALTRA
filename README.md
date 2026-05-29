# taiji-agent_ALTRA

ArchR-based single-cell pseudobulk Taiji pipeline following the [Wei Wang Lab ALTRA workflow](https://github.com/Taiji-pipeline/Taiji_ALTRA). Designed for separate-assay scRNA-seq + scATAC-seq from the same donors, using ArchR co-clustering to define pseudobulk groups and feeding raw fragment BEDs directly to Taiji (which calls peaks internally).

## Installation

```bash
git clone https://github.com/ejc043/taiji-agent_ALTRA.git
cd taiji-agent_ALTRA
claude
```

Then tell Claude:

```
Set up the ALTRA environment
```

Claude will create the `taiji-agent-altra` conda environment (~20–40 min), install ArchR (~5 min), and install SeuratData + the Azimuth PBMC reference (`pbmcref`) for celltype.l2 label transfer. Once done, activate and reopen:

```bash
micromamba activate taiji-agent-altra
claude
```

The EpiTensor HiC loop files are vendored in this repo — no download needed:

| Genome | File |
|---|---|
| hg38 | `epitensor/hg38/epitensor_loop_top10p_87311.txt` |
| hg19 | `epitensor/hg19/Top1_Perc_Epitensor.hg19.sorted.filtered.bed` |
| mm10 | `epitensor/mm10/loops.txt` |

## What you need per dataset

| File | Description |
|---|---|
| `fragments.tsv.gz` + `.tbi` | scATAC-seq fragment file (bgzipped, tabix-indexed) |
| `<sample>_labeled.h5` | scRNA-seq in HISE HDF5 format |
| Taiji binary | `bin/Taiji.1.2.0/taiji` (CentOS x86_64) |

The Azimuth PBMC reference (`pbmcref`) is downloaded automatically during `Set up the ALTRA environment` via `SeuratData::InstallData("pbmcref")` — no manual download needed.

## Running on a new dataset

Activate the env and open Claude Code:

```bash
micromamba activate taiji-agent-altra
claude
```

Then prompt Claude:

```
Run the ALTRA pseudobulk Taiji pipeline on my dataset.
  - ATAC fragments: data/<dataset>/<sample>_fragments.tsv.gz
  - RNA: data/<dataset>/<sample>_labeled.h5
  - Genome: hg38
  - Output: data/<dataset>/taiji_run/
  Log everything.
```

Claude chains three skills in order:

1. **`archr-preprocess`** — ArchR Arrow → IterativeLSI → Seurat reference transfer → addGeneIntegrationMatrix → clusKNN co-clustering → `JointCCA-cluster.rds`
2. **`pseudobulk-altra`** — one RNA TSV + one ATAC BED.gz per co-cluster → `cluster_identity.csv`
3. **`taiji-run-altra`** — builds xlsx (RNA GeneQuant + ATAC PairedEnd + HiC ChromosomeLoop), submits SLURM array → `GeneRanks.tsv` per cluster

See `skills/<skill-name>/SKILL.md` for parameter tables, output schemas, and adaptation notes.

## Runtime requirements (SLURM, hg38, ~17k ATAC cells)

| Stage | SLURM mem | CPUs | Wall time |
|---|---|---|---|
| archr-preprocess | 150G | 6 | 6h |
| pseudobulk-altra | 100G | 4 | 3h |
| taiji-run-altra (per cluster) | 30G | 4 | 24h |

## License

MIT. See `LICENSE`.

## Author

Eunice Choi (UCSD bioinformatician, ejc043@ucsd.edu).
