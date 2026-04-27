# taiji-agent

A Python package + skills bundle that wraps the [Taiji pipeline](https://github.com/Taiji-pipeline/Taiji) (Wei Wang lab, UCSD) for routine use on Mac/SLURM. The agent does not re-implement Taiji — it handles preflight validation, dataset classification, single-cell preprocessing, input construction, per-sample binary execution, and audit logging so you can hand it a directory of data and get a Taiji-ready run with minimal manual setup.

## Prerequisites

Install once on each machine where you'll run the pipeline:

| Tool | Why | Install |
|------|-----|---------|
| **git** | clone the repo | `brew install git` (Mac) / system pkg manager (Linux) |
| **python ≥ 3.10** | drives the skills + scripts | usually preinstalled; `brew install python@3.11` if not |
| **a conda solver** | builds the env from the profile-specific YAMLs | `brew install micromamba` (recommended, fastest) — `mamba` or `conda` also work |
| **Claude Code** *(optional)* | for the agentic flow described at the end | follow Anthropic's install instructions |

Linux users get a free OS-level sandbox bonus if `bubblewrap` is installed (`apt install bubblewrap` / `dnf install bubblewrap`); macOS uses `sandbox-exec` (built-in). Neither is required — soft sandboxing works without them.

## Quick start

### 1. Clone the repo

```bash
git clone https://github.com/ejc043/taiji-agent.git
cd taiji-agent
```

### 2. Auto-detect what your data needs and install

`bin/auto-install.sh` runs `detect-dataset-type` on your data directory, picks the right profile (`base` for bulk-only, `sc` for single-cell), and installs everything end-to-end:

- conda env named **`taiji-agent`** with the right packages for your dataset
- GitHub-only R packages (SeuratDisk, MuDataSeurat) if the SC profile is needed
- The Taiji binary auto-selected by `--system` (`centos`, `ubuntu`, `macos`)

```bash
# Drop your dataset under data/<your_dataset>/ first, then:
bash bin/auto-install.sh --data-dir data/ --system macos
# Linux: --system centos  or  --system ubuntu
```

What this does, in order:

1. Classifies `data/` as bulk / single-cell / mixed (detect-dataset-type).
2. Picks the profile:
   - **bulk** → `base` (~500 MB, ~5 min) — Python + xlsx + MACS3 only
   - **single-cell** → `sc` (~5 GB, ~25 min) — adds R + Seurat/Signac on top of base
   - **mixed** → `sc` (covers both branches)
3. Creates the `taiji-agent` conda env from the matching `environment.<profile>.yml`. **Skips this step if the env already has that profile installed** (idempotent — re-running on a fresh clone vs an existing env is fast).
4. Runs `bin/postinstall.R` if SC profile is in scope (installs SeuratDisk + MuDataSeurat from GitHub).
5. Detects pre-existing Taiji binaries in `binaries/` by name; if a matching one is present, just creates a symlink. Otherwise downloads the Linux binary from the Taiji GitHub release for `--system centos|ubuntu`. macOS prints manual-download instructions because the asset filename varies by macOS version.

### 3. Stage reference data (one-time per machine per genome)

```bash
micromamba activate taiji-agent
python skills/fetch-references/scripts/fetch.py \
  --genome hg38 --output dependencies_data/ --update-genomes-yml
```

Downloads to `dependencies_data/hg38/`:
- `genome.fa` (GENCODE primary assembly, ~3 GB uncompressed)
- `genes.gtf` (GENCODE v45)
- `HOCOMOCOv11_human.meme` (TF binding motifs)

Idempotent — files within ±5% of expected size are skipped. Pass `--check` to verify without downloading. Pass `--force` to re-fetch.

### 4. Verify the install

```bash
bash bin/doctor.sh --profile base       # or --profile sc / full
```

Should report PASS for every skill in the requested profile.

## Required input formats

### Bulk data

Place one file per sample per assay under a single flat directory. File-naming convention is up to you — you supply the pattern to `build-taiji-input` via `--rna-pattern`, `--atac-pattern`, and `--hic-pattern`.

| Assay | Extension | Format | Schema |
|-------|-----------|--------|--------|
| RNA-seq | `.tsv` | GeneQuant | 2-column, **no header**: `gene_symbol<TAB>expression_value` |
| ATAC-seq | `.narrowPeak` | ENCODE BED6+4 | 10-column standard MACS2/MACS3 output |
| HiC *(optional)* | `.bedpe` | ChromosomeLoop | 6-column: `chr1 start1 end1 chr2 start2 end2` |

RNA and ATAC are both required per sample; HiC is optional but improves TF→target edge accuracy.

`detect-dataset-type` will emit `MISSING REQUIRED` warnings for any absent modality when it classifies a directory as bulk. `build-taiji-input` will error (strict mode, default) or warn-and-skip (--no-strict) if a pattern resolves to a path that doesn't exist.

### Single-cell data

| Extension | Source | Loader | Modality |
|-----------|--------|--------|----------|
| `.rds` | Seurat (v5) or SingleCellExperiment | direct R `readRDS` | RNA-only / ATAC-only / multiome / separate-assay |
| `.h5ad` | AnnData (Python) | via `SeuratDisk::Convert` | typically RNA-only or separate-assay |
| `.h5mu` | MuData (Python) | via `MuDataSeurat` | always treated as multiome |

**Required companion for any ATAC processing:**

`fragments.tsv.gz` — bgzipped, Tabix-indexed scATAC fragments file (10x cellranger-arc output or Signac-compatible). MACS2/MACS3 reads barcode-filtered fragments per cluster. Without it, only RNA pseudobulks can be generated; pass `--rna-only` to `pseudobulk-construct` to suppress the error.

`detect-dataset-type` will emit a `MISSING RECOMMENDED` warning when it sees SC object files but no fragments file in the same tree.

## Running Taiji on a dataset

### Option A — direct script invocation

```bash
# Scaffold a run from the demo template:
cp -r runs/RA_OA_chr22_demo runs/<your_run_name>
rm -rf runs/<your_run_name>/{Input,Output,log,taiji_input.xlsx,taiji_config.yml,taiji_config_files.txt}

# Edit samples.csv with your sample IDs + cohorts:
$EDITOR runs/<your_run_name>/samples.csv
# Format:
#   group,cohort,rep
#   SAMPLE_01,COHORT_A,1
#   SAMPLE_02,COHORT_B,1

# If your filename pattern differs from the demo, edit bin/run-taiji.sh's
# --rna-pattern / --atac-pattern / --hic-pattern (currently set for CHEM280).

bash bin/run-taiji.sh runs/<your_run_name>
```

Outputs land at `runs/<your_run_name>/Output/Partial/<sample>/GeneRanks.tsv` (TF rankings, one per sample) and `Network/<sample>/edges_combined.csv` (TF→target edges).

### Option B — agentic flow via Claude Code

From your `taiji-agent/` directory after activating the env:

```bash
micromamba activate taiji-agent
claude
```

Then prompt Claude with something like:

```
Run Taiji on the dataset in `data/<your_dataset>/`. The samples are
<sample_id_1>, <sample_id_2>, ... — bulk RNA-seq + ATAC-seq + (optionally)
HiC, hg38. Walk through prep, classification, input construction, and the
per-sample binary execution. Validate outputs at the end. Log everything to
the workflow log so I can audit later.
```

Claude reads `CLAUDE.md` (auto-loaded from cwd), sees the six skills + the verified RA_11 baseline, and chains:

1. `detect-dataset-type` → classify
2. `fetch-references` (if not already staged) → reference data
3. `pseudobulk-construct` (only if SC) → per-cluster pseudobulks
4. `build-taiji-input` → xlsx
5. `taiji-runner` → per-sample TSV/config + per-sample Taiji invocation
6. Output validation: confirms every sample produced `GeneRanks.tsv`, summarizes top TFs

Every stage is logged to `runs/<name>/log/<run_id>.{md,jsonl}`.

## Repository layout

```
taiji-agent/
├── README.md                       this file
├── CLAUDE.md                       loaded by Claude on every session — full skills catalog + conventions
├── bin/                            install + run scripts
│   ├── auto-install.sh             agent-driven: data → profile → install
│   ├── install.sh                  profile-aware: --profile {base|sc|dev|full}
│   ├── install-taiji.sh            per-system Taiji binary downloader / detector
│   ├── postinstall.R               GitHub-only R packages (SeuratDisk, MuDataSeurat)
│   ├── run-taiji.sh                per-run executor (delegates per-sample to taiji-runner)
│   ├── sandbox-run.sh              workspace-bounded execution wrapper
│   ├── preflight-xlsx.py           verifies all paths in xlsx exist before launch
│   └── doctor.sh                   --profile filter; verifies every dep across skills
├── skills/                         six production skills (see CLAUDE.md for catalog)
│   ├── detect-dataset-type/        classifies bulk/SC/mixed; SC sub-modality
│   ├── fetch-references/           idempotent GENCODE + HOCOMOCO downloader
│   ├── build-taiji-input/          produces Taiji input xlsx
│   ├── pseudobulk-construct/       SC → bulk-Taiji bridge (Seurat WNN + per-cluster MACS3)
│   ├── taiji-runner/               per-sample Taiji 1.3.0 orchestrator
│   └── workflow-log/               per-run audit log auto-attached by every skill
├── environment.base.yml            ~500 MB Python + xlsx + MACS3
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

## Sandboxing

Every command the agent invokes runs through `bin/sandbox-run.sh`, which:

- Refuses to run from outside the workspace (`cwd` must be inside `REPO_ROOT`)
- Pins `TMPDIR` to `<REPO_ROOT>/tmp/sandbox-<pid>` so even Taiji's scratch I/O stays in-tree
- Auto-engages `bwrap` (Linux) or `sandbox-exec` (macOS) if available — kernel-enforced filesystem isolation
- Falls back to soft enforcement when no OS-level tool is present

Add `--strict-sandbox` to `bin/run-taiji.sh` (or set `TAIJI_STRICT_SANDBOX=1`) to refuse runs unless an OS-level sandbox is available. Use `--no-sandbox` only for debugging.

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
- `bin/run-taiji.sh` regenerates the xlsx + per-sample configs on every invocation; outputs from previous successful samples are preserved.

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
