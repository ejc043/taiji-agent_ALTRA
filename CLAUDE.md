# Taiji agent — repository guide for Claude

## What this repo is

A Python package + plugin-style skills bundle that wraps the [Taiji pipeline](https://github.com/Taiji-pipeline/Taiji) (Wei Wang lab, UCSD) for routine use on a SLURM HPC. The agent does not re-implement Taiji — it handles preflight validation, input construction, single-cell preprocessing, and submission so the user can hand it a directory of data and get a Taiji-ready run with minimal manual setup.

User: Eunice Chen (UCSD bioinformatician, ejc043@ucsd.edu). Practitioner-level — appreciates context, summaries, and explicit limitations. Not a beginner coder.

## Top-level layout

```
src/taiji_agent/        Python package (CLI scaffold; cli, config, detect, samplesheet, validate)
skills/                 Plugin-style skills (the bulk of the project's logic lives here)
  build-taiji-input/        produces the bulk-Taiji input xlsx
  detect-dataset-type/      classifies a directory as bulk/SC/mixed/unknown + SC sub-modality
  fetch-references/         idempotent downloader for FASTA + GTF + MEME (per genome)
  pseudobulk-construct/     bridges single-cell objects to bulk Taiji inputs
  taiji-runner/             per-sample Taiji 1.3.0 orchestrator (xlsx → per-sample TSVs/configs → taiji run per sample)
  workflow-log/             per-run audit log (md + jsonl) auto-attached by sibling skills
  taiji/                    EMPTY placeholder
  *-workspace/              fixtures + eval scaffolding for the like-named skill
  <skill>/dependencies.yml  per-skill runtime dep declaration (read by bin/doctor.sh)
bin/                    one-command install + verifier
  install.sh                top-level installer (env + R postinstall + Taiji binary)
  install-taiji.sh          per-system Taiji binary downloader
  postinstall.R             remotes::install_github() for SeuratDisk + MuDataSeurat
  doctor.sh                 verifies every dep declared in skills/*/dependencies.yml
binaries/               where install-taiji.sh drops the Taiji executable (gitignored payload)
docs/environment.md     full env-management guide (Tier 1/2/3: install, lockfile, conda-pack)
examples/demo.taiji.input.xlsx   sample of what build-taiji-input produces
templates/              SBATCH and Jinja templates for the agent's submit step (currently empty)
environment.yml         conda env spec (Python + R + Seurat/Signac + MACS3 + xlsx deps)
conda-lock.linux-64.yml lockfile for byte-reproducible installs (generate via `conda-lock`)
pyproject.toml          declares the `taiji-agent` package + `taiji-agent` CLI entry point
```

## Skills currently built

Six production skills, all chainable end-to-end. The `taiji/` directory is reserved for a future umbrella skill but is empty today.

### 1. `detect-dataset-type` — pre-flight classifier

**Purpose:** Classify a directory as `bulk` / `single-cell` / `mixed` / `unknown` from file extensions plus 10x cellranger filename patterns. For single-cell, sub-classifies into `multiome` (RNA+ATAC same cells), `separate-assay` (different cells, needs Signac `integrate_atac` co-embedding), or `sc-undetermined`. No file opens — fast and deterministic.

**Inputs:** one or more directory or file paths.

**Outputs:** a `DetectResult` dataclass with `classification`, `sc_modality`, evidence (`bulk_files`, `sc_files`), `sizes`, `warnings`, `errors`. CLI prints text or JSON; library API returns the dataclass.

**Recognized signatures:**
- Bulk: `.tsv` (RNA-seq GeneQuant), `.narrowPeak` (ATAC-seq peaks), `.bedpe` (HiC chromatin loops).
- Single-cell extensions: `.h5ad` (AnnData), `.rds` (Seurat / SingleCellExperiment), `.h5mu` (MuData — always multiome).
- Single-cell filename overrides (take priority over bulk extensions): `fragments.tsv`, `barcodes.tsv`, `features.tsv`, `genes.tsv`, `matrix.mtx`, `filtered_feature_bc_matrix.h5`, `filtered_peak_bc_matrix.h5`, `raw_feature_bc_matrix.h5`. Case-insensitive substring match; `.gz`-aware.

**Modality tier logic** (first match wins): `.h5mu` → multiome; cellranger-arc structural signature (`filtered_feature_bc_matrix.h5` + `atac_fragments.tsv`) → multiome; explicit multiome filename token (`multiome`, `multi_omic`, `multimodal`, `arc`) → multiome; both RNA hints AND ATAC hints across filenames → separate-assay; otherwise sc-undetermined.

**Used by:** `build-taiji-input` (refuses SC inputs); `pseudobulk-construct` (refuses bulk inputs and uses `sc_modality` to pick clustering signal).

**Status:** verified across 11 fixtures. Eval loop ran with 100%/100% pass-rate tie under lenient grader, 4× lower variance with-skill.

### 2. `build-taiji-input` — bulk Taiji xlsx builder

**Purpose:** Construct the input xlsx (`Active` + `active_metadata` sheets) that the bulk Taiji binary consumes. Walks a `--data-dir` for per-sample RNA `.tsv`, ATAC `.narrowPeak`, optional HiC `.bedpe`; joins against a samples manifest; emits a properly-shaped xlsx.

**Inputs:**
- `--samples` CSV/TSV with at least: `name, cohort, group, rna_seq, atac_seq, [hic, genome]`.
- `--data-dir` path to a directory containing the per-sample files.
- `--genome` tag (hg19, hg38, mm10, mm39).
- Optional `--rna-pattern` / `--atac-pattern` / `--hic-pattern` filename globs; `--no-hic`; `--strict`/`--no-strict`; `--skip-data-type-check`.

**Outputs:** an `.xlsx` at the user-specified path with two sheets matching Taiji's expected schema.

**Pre-flight gate:** soft-imports `detect-dataset-type`. If the data dir is `single-cell`, the gate refuses and surfaces the `sc_modality` in the error with modality-specific guidance (multiome → sc-Taiji; separate-assay → Signac `integrate_atac`; sc-undetermined → ask user to disambiguate). If the sibling skill isn't installed, the gate degrades to a warning.

**Status:** eval loop closed at 100% pass rate with-skill vs 11.9% without-skill — strongest delta in the project. Verified on 8 fixture configurations.

### 3. `pseudobulk-construct` — single-cell → bulk-Taiji bridge

**Purpose:** Convert a single-cell object (`.rds`, `.h5ad`, `.h5mu`) plus a fragments file into pseudobulk RNA GeneQuant TSVs and per-cluster narrowPeak files in the exact layout `build-taiji-input` consumes. Closes the SC → bulk Taiji loop.

**Inputs:**
- `--input` path to single-cell object.
- `--fragments` path to `fragments.tsv.gz` (required for ATAC peak calling).
- `--genome` tag.
- `--output-dir`.
- Optional: `--metadata-cols`, `--target-cluster-size` (default 200, band [100, 300]), `--min-cluster-cells` (default 20), `--clustering-signal` (`wnn`/`rna`/`atac`; auto-selected from `sc_modality`), `--transferred-label-col` (default `predicted.id`), `--peak-caller` (`macs2`/`macs3`), `--rna-only`/`--atac-only`, `--dry-run`, `--yes`.

**Pipeline (orchestrated by `pseudobulk.py`):**
1. `detect-dataset-type` gate — refuses bulk; for sc-undetermined refuses with Signac integrate_atac pointer; for separate-assay routes to ATAC-LSI clustering with required transferred-label check.
2. Dependency check (Rscript + MACS2 on PATH; degrades to warning under `--dry-run`).
3. `load_and_cluster.R` — loads object (extension dispatch: `.rds` direct, `.h5ad` via SeuratDisk, `.h5mu` via MuDataSeurat), runs Seurat/Signac WNN clustering (or RNA-only PCA / ATAC-only LSI), scale-aware resolution binary search (seed `r0 = 0.15·log2(N/1000)·sqrt(target_n_clusters/20)`, clamped to [0.05, 3.0], up to 8 iterations targeting mean cluster size in [100, 300]), drops <20-cell clusters and sub-groups, writes `clusters.csv`, `resolution_trace.json`, `groups_plan.json`, and per-group barcode files.
4. `aggregate_rna.R` — sums raw counts per `(cluster × metadata_col × value)` group, writes 2-column GeneQuant TSV (no header) per group.
5. Per-cluster MACS2 (driven from Python): `gunzip -c fragments.tsv.gz | awk <barcode-filter> | macs2 callpeak -t - -f BED --nomodel --shift -100 --extsize 200 -q 0.05 --call-summits`. Standard scATAC per-cluster recipe.
6. Emits `manifest.tsv` with the exact columns `build-taiji-input --samples` expects: `name, cohort, group, rna_seq, atac_seq, genome`. Zero-glue handoff.

**Outputs:**
```
output_dir/
├── clusters.csv                    per-cell barcode + cluster + metadata
├── resolution_trace.json           audit log of the resolution binary search
├── groups_plan.json                per-(cluster × meta_col × value) plan
├── per_cluster_barcodes/<name>.txt one barcode per line, used by MACS2
├── rna/<name>.tsv                  GeneQuant TSV per group
├── atac/<name>_peaks.narrowPeak    per-cluster MACS2 peaks
└── manifest.tsv                    feeds build-taiji-input directly
```

**Status:** Python orchestrator + gates verified across detect-dataset-type's 11 fixtures (bulk redirect, sc-undetermined refusal, multiome auto-WNN, separate-assay auto-ATAC, dep-check error path, --dry-run graceful degradation). R scripts authored against Seurat v5 / Signac ≥1.12 but **not yet syntax-checked or end-to-end run** because R isn't in this sandbox; needs first-run validation on Eunice's SLURM env. No fixture .rds/.h5ad/.h5mu files exist; no eval loop yet.

### 4. `workflow-log` — per-run audit trail

**Purpose:** Construct and update a per-run log (`.md` + `.jsonl`) in `<work_dir>/log/` capturing every stage of a Taiji-agent workflow. Sibling skills auto-attach via a small `active_run` pointer file — when a run is active, every call to `detect`, `pseudobulk`, or `build-taiji-input` writes a stage entry; when no run is active, the skills behave exactly as before (zero-cost soft import).

**Outputs (per run, in `<work_dir>/log/`):**
- `<run_id>.md` — human-readable narrative
- `<run_id>.jsonl` — one JSON object per stage, machine-parseable
- `index.jsonl` — one line per run (start + finalize), for cross-run grep
- `active_run` — text file holding the current run_id; deleted on `finalize`

`<run_id>` = `YYYYMMDD_HHMMSS_<4hex>` (sortable, collision-resistant).

**API (Python library + CLI):**
```python
from log import TaijiLog
log = TaijiLog.start(work_dir=".", genome="hg38",
                     reference_paths={"fasta": "...", "gtf": "...", "meme": "..."})
# sibling skills internally:
log = TaijiLog.attach_active(work_dir=".")     # returns None if no active run
if log: log.append_detect(detect_result)
log.finalize(status="success", duration_s=...)
```

```bash
python skills/workflow-log/scripts/log.py start --work-dir . --genome hg38
python skills/workflow-log/scripts/log.py status
python skills/workflow-log/scripts/log.py finalize --status success
```

**Per-stage fields captured:** classification + sample counts (detect); clustering signal + chosen resolution + full resolution_trace + n_clusters_kept/dropped + manifest_path (pseudobulk); xlsx_path + SHA-256 short hash + Active row count + dropped samples + cohorts (build_input); command + exit_code + duration + output_dir size + pagerank/edges file counts (taiji_run); plus a run header with host/user/conda env/genome/Taiji binary version/`bin/doctor.sh` snapshot.

**Soft-attach contract:** `TaijiLog.attach_active(work_dir)` returns a bound instance if `<work_dir>/log/active_run` exists, else `None`. Sibling skills always check for `None` before calling `append_*`. Logging failures inside the sibling skill are swallowed so a log bug never breaks a workflow.

**Status:** end-to-end smoke-tested in a clean working dir — start → detect → pseudobulk dry-run → finalize produces clean md + jsonl + index + the `active_run` pointer is correctly cleaned up on finalize.

### 5. `fetch-references` — per-genome reference-data downloader

**Purpose:** Idempotently download the reference data Taiji needs (FASTA + GTF + MEME motif file) from stable upstream hosts (GENCODE via EMBL-EBI mirror; HOCOMOCO v11 via autosome.org) into `<output_dir>/dependencies_data/`. Reads a single `reference_manifest.yml` declaring URLs, target paths, gunzip flags, and approximate uncompressed sizes per genome.

**Genomes available out of the box:** `hg38` (GENCODE v45), `hg19` (GENCODE v19), `mm10` (GENCODE M25), `mm39` (GENCODE M34). Add new ones by editing `scripts/reference_manifest.yml`.

**CLI:**
```bash
python scripts/fetch.py --list                                            # show genomes
python scripts/fetch.py --genome hg38 --output dependencies_data/ --check # report present/missing
python scripts/fetch.py --genome hg38 --output dependencies_data/         # download missing
python scripts/fetch.py --genome hg38 --output dependencies_data/ --update-genomes-yml
python scripts/fetch.py --genome hg38 --output dependencies_data/ --force # re-download
```

**Idempotency:** files present at the expected path with size within ±5% of the manifest's `approx_size_mb` are skipped. Downloads are atomic (`.part` + rename). `--force` re-downloads everything; `--check` reports without downloading; `--dry-run` resolves URLs without touching the network. The skill never deletes files it didn't write itself.

**Optional integrations:** `--update-genomes-yml` rewrites `skills/build-taiji-input/assets/genomes.yml` in place to point at the resolved fasta/gtf paths (preserving other genome entries verbatim). When workflow-log is active, the skill auto-attaches and writes a `fetch_references` stage entry with the per-file status, sizes, and target paths. `samtools faidx` is auto-built after FASTA download if samtools is on PATH; otherwise the skill prints a NOTE and the user can build the index manually later.

**Status:** scripts compile + run cleanly; `--list`, `--check`, `--dry-run` paths verified locally. Real downloads couldn't be exercised in this sandbox because the egress proxy blocks both EBI and autosome.org with 403 — on Eunice's Mac or SLURM node, the same commands work normally.

### 6. `taiji-runner` — per-sample Taiji 1.3.0 orchestrator

**Purpose:** Bridge between `build-taiji-input`'s xlsx (the human manifest) and Taiji's actual input format (per-sample TSV + per-sample YAML config). `taiji run --config <yml>` is single-sample-per-call and reads YAML/TSV, not xlsx — so a multi-sample run requires splitting the xlsx into N TSVs and N configs, then invoking Taiji N times. This skill does that.

**Mirrors the existing UCSD wrapper** at `Taiji/taiji_wrapper-uniqueGen.py` — same directory layout (`Input/<sample>_input/`, `Output/Partial/<sample>_output/`), same three-placeholder convention (`[insert_input_filepath_here]`, `[insert_output_directory_here]`, `[insert_genome_filepath_here]`), same `taiji_config_files.txt` index. Drop-in compatible.

**Inputs:**
- `--xlsx` path to build-taiji-input output (must have both `Active` and `active_metadata` sheets).
- `--template` per-sample formula config template (3 placeholders + optional `${REPO_ROOT}`).
- `--run-dir` where the per-sample tree lands.
- `--binary` path to Taiji executable (or `$TAIJI_BINARY`).
- Optional: `--threads N` (per-sample), `--parallel N` (samples concurrent), `--samples A,B`, `--prepare-only`, `--continue-on-error`, `--json`.

**Outputs (under `--run-dir`):**
```
Input/<sample>_input/<sample>_input.tsv     # filtered Active sheet
Input/<sample>_input/<sample>_config.yml    # template materialized for this sample
Output/Partial/<sample>_output/             # taiji's cwd; pagerank.tsv lands here
taiji_config_files.txt                      # newline list (for SLURM-array compat)
log/<sample>.taiji.{stdout,stderr}          # per-sample stream capture
```

**Output validation (defense against silent failures):** Taiji 1.3.0 sometimes exits 0 even when its workflow errored (`Program exits with errors` printed to stderr but exit code 0). The runner cross-checks: a sample is `ok` only if `exit_code == 0` AND at least one `*_pagerank.tsv` exists in its output dir. Otherwise `error: "Taiji exited 0 but produced no pagerank files — silent workflow failure."`

**Workflow-log integration:** each per-sample run becomes a separate `taiji_run.<sample>` log entry with full metadata (config_path, tsv_path, output_dir, exit_code, duration_s, n_pagerank_files, n_edges_files, stderr_tail). Soft-attach pattern; logging never breaks the runner.

**Status:** Production-validated. `--prepare-only` end-to-end verified in the sandbox (per-sample tree, configs, TSVs all generated correctly). Real `taiji run` execution verified on Eunice's Mac for sample RA_11 (4 min 44 s wall time, "Program finishes successfully"). Architectural constraint: the Taiji binary is x86_64 ELF / Mach-O; the sandbox is aarch64 with no x86_64 emulator (`Exec format error`), so Taiji **always** runs on Mac/SLURM, never in the sandbox.

**Verified Taiji 1.3.0 output schema** (per validated RA_11 run):
```
Output/Partial/<sample>_output/
├── GeneRanks.tsv               # headline output: per-TF PageRank score (765 rows for HOCOMOCO v11 human)
├── GeneRanks_PValues.tsv       # paired p-values for each TF rank
├── sciflow.db                  # workflow state DB (presence = clean checkpoint)
├── Network/<sample>/
│   ├── edges_combined.csv      # TF -> target edges, Neo4j CSV: :START_ID,:END_ID,weight,:TYPE
│   ├── edges_binding.csv       # raw binding edges (always larger than combined)
│   └── nodes.csv               # graph nodes
├── ATACSeq/
│   ├── openChromatin.bed.gz
│   └── TFBS/                   # motif-site partitions
├── RNASeq/
│   ├── expression_profile.tsv
│   └── RNA-<sample>_rep1_gene_quant.tsv
└── GENOME/genome.index         # ~3 GB; built once per genome, reused across samples
```

> **Critical naming gotcha**: Taiji 1.3.0 emits `GeneRanks.tsv`, NOT `*_pagerank.tsv`. The binary's symbol table contains stage names like `Output_Ranks` / `Output_Ranks_SC`, but those are workflow STAGES — the actual output file is `GeneRanks.tsv`. Stage names ≠ output filenames. The runner's silent-failure check looks for `GeneRanks.tsv`.

**Composition with `bin/run-taiji.sh`:** the shell wrapper at `bin/run-taiji.sh` is a thin orchestrator over this skill — it picks the binary by OS, regenerates the xlsx, runs preflight, starts/finalizes the workflow log, and delegates per-sample to `taiji-runner`. Users who want only the per-sample machinery can call `python skills/taiji-runner/scripts/run_taiji.py` directly.

## Skill composition (end-to-end pipeline)

```
    fetch-references --genome hg38 --output dependencies_data/  (one-time per site)
                              │
                              ▼
    workflow-log.start(work_dir, genome, reference_paths) ──┐
                              │                              │
                       detect-dataset-type ─┐                │
                              │             ├── auto-attach to active run
            ┌─────────────────┴───────────┐ │  via <wd>/log/active_run
            ▼                             ▼ │  pointer; emit stage entries
      classification=bulk         classification=single-cell
            │                             │
            │                    pseudobulk-construct ────────┤
            │                             │                   │
            │                         manifest.tsv            │
            │                             │                   │
            └────────►   build-taiji-input  ◄─────────────────┤
                              │                               │
                              ▼                               │
                  taiji_input.xlsx (Active + active_metadata) │
                              │                               │
                              ▼                               │
              taiji-runner: split xlsx → per-sample TSVs/configs ─┤
                              │                                  │
                              ▼                                  │
              for each sample: bulk Taiji binary               │
                  (`taiji run --config <sample>_config.yml`)   │
                              │                                  │
                              ▼                                  │
              workflow-log.finalize(status) ─────────────────────┘
```

Each step gates on the previous one and produces inputs in the next step's expected format with **no glue code in between**. The `manifest.tsv` from `pseudobulk-construct` is a drop-in for `build-taiji-input --samples`.

## External dependencies (not vendored)

Trimmed to reflect that all inputs are pre-aligned/pre-quantified — no BWA, STAR, RSEM, samtools, Picard, or Singularity needed. Per-skill dep manifests live at `skills/*/dependencies.yml`; `bash bin/doctor.sh` verifies them.

**Always required:**
- Python ≥ 3.10
- Pinned Python deps in `environment.yml`: click, pydantic, pyyaml, jinja2, rich, pandas, openpyxl
- Rscript ≥ 4.2 (only for `pseudobulk-construct`)
- R packages (bioconda): Seurat ≥ 5.0, Signac ≥ 1.12, Matrix, GenomicRanges, dplyr, optparse, jsonlite
- R packages (GitHub-only, installed by `bin/postinstall.R`): SeuratDisk (for .h5ad), MuDataSeurat (for .h5mu)
- MACS3 ≥ 3.0 (or MACS2 — `pseudobulk-construct` auto-detects whichever is on PATH)
- gunzip + awk (standard Unix)

**Taiji binary** — installed by `bin/install-taiji.sh --system <centos|ubuntu|macos>`:
- CentOS / RHEL: `taiji-CentOS-x86_64` (auto-fetched)
- Ubuntu / Debian: `taiji-Ubuntu-x86_64` (auto-fetched)
- macOS: `taiji-macOS-XX-XX` (suffix varies by macOS version; install-taiji.sh prints instructions)
- Singularity: explicitly out of scope (no official .sif published, can't install on this site's SLURM)

For Eunice's UCSD setup, the binary additionally lives at `/stg3/data1/eunice/Projects/nextflow/nf_taiji/binaries/` and may be symlinked from there if preferred.

## Environment management

Single command to set up everything:

```bash
bash bin/install.sh
micromamba activate taiji-agent
bash bin/doctor.sh           # verify every dep across all skills
```

Three reproducibility tiers (full details in `docs/environment.md`):

1. **Tier 1 — micromamba + environment.yml.** One command, ~10 min on first run. Versions pinned via `>=` constraints; transitive resolution at install time.
2. **Tier 2 — conda-lock.** Generate `conda-lock.linux-64.yml` once, commit it, then `bash bin/install.sh --use-lockfile` produces byte-identical envs across machines.
3. **Tier 3 — conda-pack tarball.** After Tier 2 works on one machine, pack to a portable `.tar.gz`. Restore on any compatible Linux node in <2 min, no network. Closest equivalent to a Docker image without needing a container runtime.

Per-skill `dependencies.yml` files declare runtime needs (binaries with version pins + alternatives, Python packages with `import_name` overrides, R packages with bioconda/github source tags, optional data files). `bin/doctor.sh --json` reads them for machine-readable status.

## Design conventions across all skills

These are load-bearing — preserve them when adding new skills or extending existing ones.

- **CLI is a thin wrapper over a pure library function.** Every skill exposes a callable Python (or R) function returning a dataclass; the CLI is just argparse + `print(result.summary())`.
- **No file content sniffing.** All classification and routing happens on extensions and filenames. The user's lab conventions are unambiguous and the cost of opening files (h5py, R runtime) is not worth the marginal robustness.
- **Strict-on-mixed by default.** Silent SC contamination of a bulk pipeline produces garbage at Taiji-invocation time. Loud-fail-early > debug-at-analysis-time.
- **Soft-import siblings.** Cross-skill dependencies (e.g. `build-taiji-input` calls `detect-dataset-type`) use `sys.path.insert(0, sibling/scripts)` and a `try: from detect import ...` so a standalone install of one skill still works without the other.
- **Reference docs in `references/`, not in `SKILL.md`.** SKILL.md stays short; deep edge-case reasoning lives in `references/<topic>.md`.
- **Every recognized format/pattern ships with a fixture and an assertion.** `<skill>-workspace/fixtures/` mirrors `evals/evals.json`. When a verification pass finds a bug, add the failing case as a fixture.
- **Always cite limitations.** Eunice prefers explicit "this heuristic can fail when X" notes over silent confidence.
- **Output formats designed for downstream chaining.** When skill A's output feeds skill B, A emits exactly what B expects with no transformation step in between.

## Memory system

Persistent memory at `/sessions/quirky-festive-mayer/mnt/.auto-memory/`. Key files:

- `MEMORY.md` — index, loaded automatically.
- `user_role.md` — Eunice's role and preferences.
- `project_taiji_agent.md` — high-level project decisions.
- `project_detect_dataset_type.md` — skill-specific design notes.
- `project_pseudobulk_construct.md` — skill-specific design notes.

Read these before making non-trivial decisions about extending the project.

## Architectural constraints (load-bearing)

- **The Taiji binary cannot be executed from this sandbox.** Sandbox is aarch64; binaries are x86_64 ELF (Linux variants) or Mach-O (macOS). No qemu-user / box64 / fex-emu emulator is installed, and `sudo apt install` is blocked by `no_new_privileges`. The kernel rejects with `Exec format error`. Practical consequence: every "run Taiji" step is either Mac-side (Eunice's machine) or SLURM-side. The agent's role from the sandbox is strictly **prepare inputs + validate outputs**.
- **Workspace folder is shared between sandbox and Mac.** Files dropped into `<workspace>/data/`, `<workspace>/dependencies_data/`, `<workspace>/runs/<name>/Output/` from either side become immediately visible to the other. Use this for handoffs: Mac runs Taiji → outputs land in workspace → sandbox can validate / summarize / log without re-running.
- **R / MACS3 are also not in the sandbox.** Only `pseudobulk-construct` needs them; the SC branch is also Mac-side.

## Not yet built (scoped but pending)

From `project_taiji_agent.md`:
- The actual Taiji CLI plugin wrapping these skills (the `skills/taiji/` directory is reserved for this).
- A SLURM submission skill that takes a `taiji_input.xlsx` + a config and emits an sbatch script.
- A cron watchdog that monitors upstream output directories (e.g. cellranger output dirs) and triggers the pipeline when new datasets land.

For `pseudobulk-construct` specifically:
- Real fixture `.rds` / `.h5ad` / `.h5mu` files.
- An eval loop matching `build-taiji-input`'s and `detect-dataset-type`'s.
- First-run validation on the SLURM env (R syntax check + end-to-end on a public PBMC multiome).

## Where to find what

| Question                                  | Where                                                                 |
|-------------------------------------------|-----------------------------------------------------------------------|
| "What does skill X do?"                   | `skills/X/SKILL.md`                                                   |
| "Why does skill X make decision Y?"       | `skills/X/references/`                                                |
| "What design decisions were locked in?"   | `/sessions/quirky-festive-mayer/mnt/.auto-memory/project_*.md`        |
| "What fixtures cover what cases?"         | `skills/X-workspace/fixtures/` + `skills/X/evals/evals.json`          |
| "How is the Python package wired up?"     | `pyproject.toml` + `src/taiji_agent/`                                 |
| "What conda env is expected?"             | `environment.yml`                                                     |
| "What does a Taiji input xlsx look like?" | `examples/demo.taiji.input.xlsx`                                      |
