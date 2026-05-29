# taiji-agent

A Python package + skills bundle that wraps the [Taiji pipeline](https://github.com/Taiji-pipeline/Taiji) (Wei Wang lab, UCSD) for routine use on Mac/SLURM. The agent does not re-implement Taiji — it handles preflight validation, dataset classification, single-cell preprocessing, input construction, per-sample binary execution, and audit logging so you can hand it a directory of data and get a Taiji-ready run with minimal manual setup.

## Prerequisites

Install once on each machine where you'll run the pipeline:

| Tool | Why | Install |
|------|-----|---------|
| **git** | clone the repo | `brew install git` (Mac) / system pkg manager (Linux) |
| **a conda solver** | builds the env from the profile-specific YAMLs | `brew install micromamba` (recommended, fastest) |
| **Claude Code** | drives the agentic pipeline | follow Anthropic's install instructions |

Python ≥ 3.10 is installed automatically inside the `taiji-agent` conda env; no system Python version requirement.

## Quick start

### 1. Clone the repo

```bash
git clone https://github.com/ejc043/taiji-agent.git
cd taiji-agent
```

### 2. Set up the environment

**Prompt-driven (recommended):** activate Claude Code, then describe your setup:

```
micromamba activate taiji-agent   # if env exists; skip on first run
claude
```

Then tell Claude:

```
Set up taiji-agent for hg38, my data is in data/my_dataset/ on centos
```

Claude will invoke the `/setup` skill, which:
1. Runs `bin/preflight.py` to check what's already installed (env profile, genome refs)
2. Classifies `data/my_dataset/` to determine whether `base` or `sc` profile is needed
3. Runs `bin/install.sh` and/or `fetch-references` only for what's missing
4. Verifies with `bin/doctor.sh`

**Direct CLI alternative** (if you know your dataset type):

```bash
bash bin/auto-install.sh --data-dir data/ --system macos --fetch-references --genome hg38
# Linux: --system centos  or  --system ubuntu
# Replace hg38 with hg19 or mm10 as needed.
```

What `auto-install.sh` does, in order:

1. Classifies `data/` as bulk / single-cell
2. Picks the profile:
   - **bulk** → `base` (~500 MB, ~5 min) — Python + input construction + MACS3 only
   - **single-cell** → `sc` (~5 GB, ~25 min) — adds R + Seurat/Signac on top of base
3. Creates the `taiji-agent` conda env from the matching `environment.<profile>.yml`. Idempotent — re-running on an existing env skips already-installed profiles.
4. Runs `bin/postinstall.R` if SC profile is in scope.
5. Downloads the Taiji binary for `--system centos|ubuntu|macos`.

### 3. Verify the install

```bash
bash bin/doctor.sh --profile base       # or --profile sc
```

Should report PASS for every skill in the requested profile.

## Required input formats

### Bulk data

Place one file per sample per assay under a single flat directory. File-naming convention is up to you — you supply the pattern to `build-taiji-input` via `--rna-pattern`, `--atac-pattern`, and `--hic-pattern`.

| Assay | Extension | Format | Schema |
|-------|-----------|--------|--------|
| RNA-seq | `.tsv` | GeneQuant | 2-column, **no header**: `gene_symbol<TAB>expression_value` |
| ATAC-seq | `.narrowPeak` | ENCODE BED6+4 | 10-column standard MACS2/MACS3 output |
| HiC *(optional)* | `.bedpe` | ChromosomeLoop | 6-column: `chr1 start1 end1 chr2 start2 end2`. If omitted, vendored EpiTensor HiC predictions for the genome are used automatically. |

RNA and ATAC are both required per sample; HiC is optional but improves TF→target edge accuracy.

### Single-cell data

| Extension | Source | Modality |
|-----------|--------|----------|
| `.rds` | Seurat (v5) | RNA-only / ATAC-only / multiome / separate-assay |
| `.h5ad` | AnnData (Python) | typically RNA-only or separate-assay |

**Required companion for any ATAC processing:** `fragments.tsv.gz` — bgzipped, Tabix-indexed scATAC fragments file (10x cellranger-arc output or Signac-compatible). Without it, only RNA pseudobulks can be generated; pass `--rna-only` to `pseudobulk-construct` to suppress the error.

HiC is not a direct input for single-cell workflows — it enters via `build-taiji-input` after pseudobulking. Vendored EpiTensor predictions are used automatically if none is supplied.

## Running Taiji on a dataset

From your `taiji-agent/` directory after activating the env:

```bash
micromamba activate taiji-agent
claude
```

Then prompt Claude with something like:

**Bulk RNA + ATAC (+ optionally HiC):**

```
Run Taiji on the dataset in `data/<your_dataset>/`, hg38.
Log everything to the workflow log so I can audit later.
```

**Single-cell, multiome (cellranger-arc — same cells in both modalities):**

```
Run Taiji on the multiome dataset at `data/<your_dataset>/`, hg38. Stratify
pseudobulks by <metadata cols, e.g. condition,donor>.
Log everything to the workflow log.
```

**Single-cell, separate-assay (one .rds RNA + one .rds ATAC, different cells):**

```
Run Taiji on the dataset in `data/<your_dataset>/`, hg38.
Log everything to the workflow log.
```

Claude reads `CLAUDE.md` (auto-loaded from cwd) and chains the eight production skills:

1. `detect-dataset-type` → classify bulk vs. SC and sub-modality
2. `fetch-references` (if not already staged) → GENCODE FASTA + GTF + CIS-BP MEME
3. `sc-qc` (SC only) → per-cell QC filters before any downstream step
4. `coembed-construct` (separate-assay SC only) → joint PCA/UMAP + de novo clusters
5. `pseudobulk-construct` (SC only) → per-cluster pseudobulks + `qc/umap.png`
6. `build-taiji-input` → Taiji input XLSX (`Active` + `active_metadata` sheets)
7. `taiji-runner` → per-sample TSV/config + per-sample Taiji binary invocation
8. Output validation: confirms every sample produced `GeneRanks.tsv`, summarizes top TFs

Every stage is logged to `runs/<name>/log/<run_id>.{md,jsonl}`.

## Repository layout

```
taiji-agent/
├── README.md                       this file
├── CLAUDE.md                       loaded by Claude on every session — full skills catalog + conventions
├── bin/                            install + run scripts
│   ├── auto-install.sh             agent-driven: data → profile → install (direct CLI)
│   ├── install.sh                  profile-aware: --profile {base|sc|dev}
│   ├── install-taiji.sh            per-system Taiji binary downloader / detector
│   ├── postinstall.R               placeholder for future GitHub-only R packages
│   ├── preflight.py                check env profiles + genome refs + data classification
│   ├── preflight-xlsx.py           verifies all paths in Taiji input XLSX before launch
│   ├── run-taiji.sh                per-run executor (delegates per-sample to taiji-runner)
│   ├── sandbox-run.sh              workspace-bounded execution wrapper
│   └── doctor.sh                   --profile filter; verifies every dep across skills
├── skills/                         eight production skills (see CLAUDE.md for catalog)
│   ├── sc-qc/                      cell-level QC filtering (RNA: nFeature/percent.mt; ATAC: TSS/nucleosome/blacklist)
│   ├── detect-dataset-type/        classifies bulk/SC; SC sub-modality
│   ├── fetch-references/           idempotent GENCODE downloader + CIS-BP staging
│   ├── coembed-construct/          separate-assay scRNA + scATAC → shared latent space
│   ├── build-taiji-input/          produces Taiji input XLSX
│   ├── pseudobulk-construct/       SC → bulk-Taiji bridge (Seurat WNN + per-cluster MACS3)
│   ├── taiji-runner/               per-sample Taiji 1.3.0 orchestrator
│   └── workflow-log/               per-run audit log auto-attached by every skill
├── environment.base.yml            ~500 MB Python + input construction + MACS3
├── environment.sc.yml              +3-4 GB R + Seurat + Signac (additive on base)
├── environment.dev.yml             +500 MB pytest + ruff + mypy
├── docs/
│   ├── environment.md              full env-management guide (lockfile, conda-pack)
│   └── workflow_diagram.md         visual pipeline diagram
└── runs/                           per-dataset outputs (gitignored)
    └── <name>/
        ├── code/                   you author everything here
        │   ├── samples.csv             per-sample manifest
        │   └── taiji_config.template.yml   per-project Taiji config (3 path placeholders)
        ├── Input/<sample>_input/   per-sample TSV + materialized config (generated)
        ├── Output/Partial/<sample>_output/   GeneRanks.tsv, Network/, etc. (generated)
        └── log/                    audit trail (md + jsonl per run; generated)
```

## Single-cell pseudobulk Taiji (ALTRA protocol)

A second, ArchR-based pipeline for scRNA-seq + scATAC-seq data that follows the Wei Wang Lab's published ALTRA workflow (`Taiji_ALTRA/scripts/prepare_taiji_input.r`). Use this instead of the standard skill chain when:

- You have **separate-assay** scRNA-seq + scATAC-seq from the **same donors** (not multiome same-cell)
- You have a **barcode quality whitelist** from an upstream QC pipeline (e.g. population demultiplexing with PassQC + singlet calls)
- You want **ArchR-based co-clustering** (rather than Seurat WNN) to define pseudobulk groups

### Additional setup (beyond the standard taiji-agent env)

Create the `taiji-agent-altra` conda environment — this installs R, Seurat, all Bioconductor dependencies, `macs2`, `bedGraphToBigWig`, and the Python helper packages in one step. Then run the post-install script to pull ArchR from GitHub (the only dependency not on conda):

```bash
# Step 1 — create env (~20–40 min, one-time)
micromamba create -f environment.altra.yml

# Step 2 — install ArchR itself from GitHub (~5 min, one-time)
micromamba activate taiji-agent-altra
Rscript bin/postinstall_altra.R
```

No Singularity container required.

The EpiTensor HiC loop files are **vendored in this repo** — no download needed:

| Genome | File |
|---|---|
| hg38 | `epitensor/hg38/epitensor_loop_top10p_87311.txt` |
| hg19 | `epitensor/hg19/Top1_Perc_Epitensor.hg19.sorted.filtered.bed` |
| mm10 | `epitensor/mm10/loops.txt` |

### What you need per dataset

| File | Description |
|---|---|
| `fragments.tsv.gz` + `.tbi` | scATAC-seq fragment file (bgzipped, tabix-indexed) |
| `<sample>_labeled.h5` | scRNA-seq in HISE HDF5 format |
| `pbmc_reference.rds` | Seurat PBMC reference (converted from `pbmc_multimodal.h5seurat`, Platelet-filtered) |
| `raw_individual_atac_metadata_backup.rds` | ATAC barcode whitelist with `PassQC` + `singlet` columns |
| `raw_individual_rna_metadata_backup.rds` | RNA barcode metadata |
| Singularity container | `archr_1.0.3.sif` (ArchR 1.0.3, Seurat 5.3.0) |
| Taiji binary | `bin/Taiji.1.2.0/taiji` (CentOS-compatible) |

### How to prompt Claude

```
Run the ALTRA pseudobulk Taiji pipeline on my dataset.
  - ATAC fragments: data/<dataset>/<sample>_fragments.tsv.gz
  - RNA: data/<dataset>/<sample>_labeled.h5
  - PBMC reference: data/<dataset>/pbmc_reference.rds
  - ATAC whitelist: backup_meta/raw_individual_atac_metadata_backup.rds (sample id: <your-atac-id>)
  - RNA whitelist: backup_meta/raw_individual_rna_metadata_backup.rds (sample id: <your-rna-id>)
  - Genome: hg38
  - Container: /path/to/archr_1.0.3.sif
  - Output: data/<dataset>/taiji_run/
  Log everything.
```

Claude will chain three skills in order:

1. **`archr-preprocess`** — ArchR Arrow → IterativeLSI → Seurat reference transfer → addGeneIntegrationMatrix → clusKNN co-clustering → `JointCCA-cluster.rds`
2. **`pseudobulk-altra`** — one RNA TSV + one ATAC BED.gz per co-cluster → `cluster_identity.csv`
3. **`taiji-run-altra`** — builds xlsx (RNA GeneQuant + ATAC PairedEnd + HiC ChromosomeLoop), submits SLURM array → `GeneRanks.tsv` per cluster

See `skills/<skill-name>/SKILL.md` for parameter tables, output schemas, and adaptation notes.

### Key differences from standard skill chain

| Standard pipeline | ALTRA pipeline |
|---|---|
| Seurat WNN clustering | ArchR IterativeLSI + Seurat CCA co-clustering |
| Pre-called MACS3 peaks (narrowPeak) | Raw fragment BED.gz — Taiji calls peaks internally |
| `pseudobulk-construct` + `build-taiji-input` | `pseudobulk-altra` + `taiji-run-altra` |
| HiC optional (auto-vendored) | HiC required (EpiTensor loops file) |
| Works on Mac or SLURM | Requires SLURM + Singularity container |

### Runtime requirements (SLURM, hg38, ~17k ATAC cells)

| Stage | SLURM mem | CPUs | Wall time |
|---|---|---|---|
| archr-preprocess | 150G | 6 | 6h |
| pseudobulk-altra | 100G | 4 | 3h |
| taiji-run-altra (per cluster) | 30G | 4 | 24h |

macs2 and bedGraphToBigWig must be in PATH for the Taiji SLURM tasks. See `skills/taiji-run-altra/SKILL.md` for paths.

## Common workflows

```bash
# Re-run a specific failed sample
bash bin/run-taiji.sh runs/<name> --samples FAILED_SAMPLE

# Run multiple samples concurrently (each uses --threads cores)
PARALLEL=2 THREADS=4 bash bin/run-taiji.sh runs/<name>

# Just prepare per-sample inputs without invoking Taiji
bash bin/run-taiji.sh runs/<name> --prepare-only

# Check what's installed vs. what's needed (without installing anything)
python bin/preflight.py --data-dir data/<dir> --genome hg38

# Verify everything is wired up
bash bin/doctor.sh --profile sc

# Inspect a previous run
cat runs/<name>/log/<run_id>.md
ls runs/<name>/Output/Partial/<sample>/GeneRanks.tsv
```

## Idempotency

All install and fetch steps are designed to be safely re-runnable:

- `bin/install.sh` tracks installed profiles in `<env>/.taiji-agent-profiles` with a hash of the env-file. Re-runs on an unchanged env are <1 s no-ops.
- `bin/install-taiji.sh` detects pre-existing binaries by name and just (re)creates the symlink — no re-download.
- `bin/preflight.py` reports current state without taking any action; use it to inspect before installing.
- `fetch-references --check` reports what's present without touching the network. Real fetches skip files already at the expected path with the right size.
- `bin/run-taiji.sh` regenerates the input XLSX + per-sample configs on every invocation; outputs from previous successful samples are preserved.

## Documentation

- **CLAUDE.md** — canonical skills catalog + conventions, auto-loaded by Claude. Read this first if you're extending the agent.
- **docs/environment.md** — full env-management guide: profiles, lockfiles, conda-pack tarballs.
- **docs/workflow_diagram.md** — visual pipeline diagram.
- **skills/<name>/SKILL.md** — per-skill documentation (triggers, I/O, conventions).
- **skills/<name>/references/** — deep dives on edge cases (clustering strategy, peak calling, output format, etc.).

## Architectural constraints

- **The Taiji binary is x86_64** (Linux ELF / macOS Mach-O). It cannot run on aarch64 without an emulator. Use the macOS binary on Apple Silicon under Rosetta if needed.
- **All Taiji invocations happen on your machine** (Mac or SLURM). The agent's role is prep + validation; it does not re-implement Taiji's compute.
- **Reference data is not in git** (~5 GB per genome). Always re-download via `fetch-references` per machine.
- **`runs/*/Output/`, `runs/*/Input/`, `runs/*/log/` are gitignored.** Per-run artifacts shouldn't pollute the repo history.

## License

MIT. See `LICENSE`.

## Author

Eunice Choi (UCSD bioinformatician, ejc043@ucsd.edu).
