# taiji-agent

A Python package + skills bundle that wraps the [Taiji pipeline](https://github.com/Taiji-pipeline/Taiji) (Wei Wang lab, UCSD) for routine use on Mac/SLURM. The agent does not re-implement Taiji — it handles preflight validation, dataset classification, single-cell preprocessing, input construction, per-sample binary execution, and audit logging so you can hand it a directory of data and get a Taiji-ready run with minimal manual setup.

## Prerequisites

Install once on each machine where you'll run the pipeline:

| Tool | Why | Install |
|------|-----|---------|
| **git** | clone the repo | `brew install git` (Mac) / system pkg manager (Linux) |
| **python ≥ 3.10** | drives the skills + scripts | usually preinstalled; `brew install python@3.11` if not |
| **a conda solver** | builds the env from the profile-specific YAMLs | `brew install micromamba` (recommended, fastest) |
| **Claude Code** | drives the agentic pipeline | follow Anthropic's install instructions |


## Quick start

### 1. Clone the repo

```bash
git clone https://github.com/ejc043/taiji-agent.git
cd taiji-agent
```

### 2. Auto-detect what your data needs and install

`bin/auto-install.sh` runs `detect-dataset-type` on your data directory, picks the right profile (`base` for bulk-only, `sc` for single-cell), and installs everything end-to-end:

- conda env named **`taiji-agent`** with the right packages for your dataset
- R packages needed for the SC profile
- The Taiji binary auto-selected by `--system` (`centos`, `ubuntu`, `macos`)

```bash
# Drop your dataset under data/<your_dataset>/ first, then:
bash bin/auto-install.sh --data-dir data/ --system macos --fetch-references --genome hg38
# Linux: --system centos  or  --system ubuntu
# Replace hg38 with hg19 or mm10 as needed.
# Omit --fetch-references to skip reference staging (re-run install.sh with the flag later).
```

What this does, in order:

1. Classifies `data/` as bulk / single-cell 
2. Picks the profile:
   - **bulk** → `base` (~500 MB, ~5 min) — Python + parent input + MACS3 only
   - **single-cell** → `sc` (~5 GB, ~25 min) — adds R + Seurat/Signac on top of base
3. Creates the `taiji-agent` conda env from the matching `environment.<profile>.yml`. **Skips this step if the env already has that profile installed** (idempotent — re-running on a fresh clone vs an existing env is fast).
4. Runs `bin/postinstall.R` if SC profile is in scope (no-op currently; kept for future GitHub-only R packages).
5. Detects pre-existing Taiji binaries in `binaries/` by name; if a matching one is present, just creates a symlink. Otherwise downloads the Linux binary from the Taiji GitHub release for `--system centos|ubuntu`. macOS prints manual-download instructions because the asset filename varies by macOS version.

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
| HiC *(optional)* | `.bedpe` | ChromosomeLoop | 6-column: `chr1 start1 end1 chr2 start2 end2`. If omitted, vendored EpiTensor HiC predictions for the genome (hg19, hg38, mm10) are used automatically. |

RNA and ATAC are both required per sample; HiC is optional but improves TF→target edge accuracy.

`detect-dataset-type` will emit `MISSING REQUIRED` warnings for any absent modality when it classifies a directory as bulk. `build-taiji-input` will error (strict mode, default) or warn-and-skip (--no-strict) if a pattern resolves to a path that doesn't exist.

### Single-cell data

| Extension | Source | Loader | Modality |
|-----------|--------|--------|----------|
| `.rds` | Seurat (v5) or SingleCellExperiment | direct R `readRDS` | RNA-only / ATAC-only / multiome / separate-assay |
| `.h5ad` | AnnData (Python) | via `SeuratDisk::Convert` | typically RNA-only or separate-assay |

HiC is not a direct input for single-cell workflows — it enters via `build-taiji-input` after pseudobulking. If no HiC is supplied at that stage, vendored EpiTensor predictions for the genome are used automatically (same fallback as bulk).

**Required companion for any ATAC processing:**

`fragments.tsv.gz` — bgzipped, Tabix-indexed scATAC fragments file (10x cellranger-arc output or Signac-compatible). MACS2/MACS3 reads barcode-filtered fragments per cluster. Without it, only RNA pseudobulks can be generated; pass `--rna-only` to `pseudobulk-construct` to suppress the error.

`detect-dataset-type` will emit a `MISSING RECOMMENDED` warning when it sees SC object files but no fragments file in the same tree.

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
Log everything to
the workflow log so I can audit later.
```

**Single-cell, multiome (cellranger-arc — same cells in both modalities):**

```
Run Taiji on the multiome dataset at `data/<your_dataset>/`, hg38. Stratify
pseudobulks by <metadata cols, e.g. condition,donor>.
Log everything to the workflow log.
```

**Single-cell, separate-assay (one .rds RNA + one .rds ATAC, different cells):**

```
Run Taiji in `data/<your_dataset>/`, hg38. 
Log everything to the workflow log.
```

Claude reads `CLAUDE.md` (auto-loaded from cwd), sees the seven skills + the verified RA_11 baseline, and chains:

1. `detect-dataset-type` → classify
2. `fetch-references` (if not already staged) → reference data
3. `coembed-construct` (only for separate-assay SC) → coembed.rds with shared PCA/UMAP and de novo clusters
4. `pseudobulk-construct` (only if SC) → per-cluster pseudobulks + `qc/umap.png`
5. `build-taiji-input` → parent input
6. `taiji-runner` → per-sample TSV/config + per-sample Taiji invocation
7. Output validation: confirms every sample produced `GeneRanks.tsv`, summarizes top TFs

Every stage is logged to `runs/<name>/log/<run_id>.{md,jsonl}`.

## Repository layout

```
taiji-agent/
├── README.md                       this file
├── CLAUDE.md                       loaded by Claude on every session — full skills catalog + conventions
├── bin/                            install + run scripts
│   ├── auto-install.sh             agent-driven: data → profile → install
│   ├── install.sh                  profile-aware: --profile {base|sc|dev}
│   ├── install-taiji.sh            per-system Taiji binary downloader / detector
│   ├── postinstall.R               placeholder for future GitHub-only R packages
│   ├── run-taiji.sh                per-run executor (delegates per-sample to taiji-runner)
│   ├── sandbox-run.sh              workspace-bounded execution wrapper
│   ├── preflight-xlsx.py           verifies all paths in parent input exist before launch
│   └── doctor.sh                   --profile filter; verifies every dep across skills
├── skills/                         six production skills (see CLAUDE.md for catalog)
│   ├── detect-dataset-type/        classifies bulk/SC; SC sub-modality
│   ├── fetch-references/           idempotent GENCODE downloader + CIS-BP staging
│   ├── build-taiji-input/          produces Taiji input parent input
│   ├── pseudobulk-construct/       SC → bulk-Taiji bridge (Seurat WNN + per-cluster MACS3)
│   ├── taiji-runner/               per-sample Taiji 1.3.0 orchestrator
│   └── workflow-log/               per-run audit log auto-attached by every skill
├── environment.base.yml            ~500 MB Python + parent input + MACS3
├── environment.sc.yml              +3-4 GB R + Seurat + Signac (additive)
├── environment.dev.yml             +500 MB pytest + ruff + mypy
├── docs/
│   ├── environment.md              full env-management guide (lockfile, conda-pack)
│   ├── workflow_diagram.md         visual pipeline diagram
│   └── workflow_slide.txt          16:9 presentation slide
└── runs/                           per-dataset outputs (gitignored)
    └── <name>/
        ├── samples.csv             you author this
        ├── taiji_config.template.yml   per-project formula (3 placeholders)
        ├── Input/<sample>_input/   per-sample TSV + materialized config
        ├── Output/Partial/<sample>_output/   GeneRanks.tsv, Network/, etc.
        └── log/                    audit trail (md + jsonl per run)
```

## Common workflows

```bash
# Re-run a specific failed sample
bash bin/run-taiji.sh runs/<name> --samples FAILED_SAMPLE

# Run multiple samples concurrently (each uses --threads cores)
PARALLEL=2 THREADS=4 bash bin/run-taiji.sh runs/<name>

# Just prepare per-sample inputs without invoking Taiji
bash bin/run-taiji.sh runs/<name> --prepare-only

# Verify everything is wired up
bash bin/doctor.sh --profile sc

# Inspect a previous run
cat runs/<name>/log/<run_id>.md
ls runs/<name>/Output/Partial/<sample>/GeneRanks.tsv
```

## Idempotency

All install steps are designed to be safely re-runnable:

- `bin/install.sh` tracks installed profiles in `<env>/.taiji-agent-profiles` with hash of the env-file. Re-runs on an unchanged env are <1 s no-ops.
- `bin/install-taiji.sh` detects pre-existing binaries by name and just (re)creates the symlink — no re-download.
- `fetch-references --check` reports what's present without touching the network. Real fetches skip files already at the expected path with the right size.
- `bin/run-taiji.sh` regenerates the parent input + per-sample configs on every invocation; outputs from previous successful samples are preserved.

## Documentation

- **CLAUDE.md** — the canonical skills catalog + conventions, auto-loaded by Claude. Read this first if you're extending the agent.
- **docs/environment.md** — full env-management guide: profiles, lockfiles, conda-pack tarballs.
- **docs/workflow_diagram.md** — visual pipeline diagram.
- **docs/workflow_slide.txt** — 16:9 presentation slide for talks.
- **skills/<name>/SKILL.md** — per-skill documentation (triggers, I/O, conventions).
- **skills/<name>/references/** — deep dives on edge cases (clustering strategy, peak calling, output format, etc.).

## Architectural constraints

- **The Taiji binary is x86_64** (Linux ELF / macOS Mach-O). It cannot run on aarch64 systems without an emulator. Use the macOS binary on Apple Silicon under Rosetta if needed.
- **All Taiji invocations happen on your machine** (Mac or SLURM). The agent's role is prep + validation; it does not re-implement Taiji's compute.
- **Reference data is not in git** (~5 GB per genome). Always re-download via `fetch-references` per machine.
- **`runs/*/Output/`, `runs/*/Input/`, `runs/*/log/` are gitignored.** Per-run artifacts shouldn't pollute the repo history.

## License

MIT. See `LICENSE`.

## Author

Eunice Chen (UCSD bioinformatician, ejc043@ucsd.edu).
