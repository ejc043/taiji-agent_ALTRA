# Taiji agent — repository guide for Claude

> **Operational QC and skill invocation guidance lives in `skills/CLAUDE.md`.** Read it before running any single-cell skill. It covers mandatory sc-qc checks, chooseR interpretation, and common failure modes found in live runs.

## What this repo is

A Python package + plugin-style skills bundle that wraps the [Taiji pipeline](https://github.com/Taiji-pipeline/Taiji) (Wei Wang lab, UCSD) for routine use on a SLURM HPC and on a Mac. The agent does not re-implement Taiji — it handles preflight validation, input construction, single-cell preprocessing, and submission so the user can hand it a directory of data and get a Taiji-ready run with minimal manual setup.

User: practitioner-level bioinformatician. Appreciates context, summaries, and explicit limitations. Not a beginner coder.

## Top-level layout

```
src/taiji_agent/        Python package (CLI scaffold; cli, config, detect, samplesheet, validate)
skills/                 Plugin-style skills (the bulk of the project's logic lives here)
  build-taiji-input/        produces the bulk-Taiji input xlsx
  coembed-construct/        co-embeds separate-assay scRNA + scATAC into shared latent space
  detect-dataset-type/      classifies a directory as bulk/SC/mixed/unknown + SC sub-modality
  fetch-references/         idempotent stager for FASTA + GTF + MEME (per genome)
  pseudobulk-construct/     bridges single-cell objects to bulk-Taiji inputs
  sc-qc/                    cell-level QC filtering (nFeature/percent.mt for RNA; TSS/nucleosome/blacklist for ATAC)
  taiji-runner/             per-sample Taiji 1.3.0 orchestrator (xlsx → per-sample TSVs/configs → taiji run per sample)
  workflow-log/             per-run audit log (md + jsonl) auto-attached by sibling skills
  taiji/                    EMPTY placeholder for future umbrella skill
  *-workspace/              fixtures + eval scaffolding for the like-named skill
  <skill>/dependencies.yml  per-skill runtime dep declaration (read by bin/doctor.sh)
bin/                    one-command install + verifier
  install.sh                top-level installer (env + R postinstall + Taiji binary)
  install-taiji.sh          per-system Taiji binary downloader
  postinstall.R             placeholder for future GitHub-only R packages
  doctor.sh                 verifies every dep declared in skills/*/dependencies.yml
  run-taiji.sh              top-level Taiji executor (binary picker + preflight + workflow-log + delegate)
  sandbox-run.sh            workspace-bounded execution wrapper
binaries/               where install-taiji.sh drops the Taiji executable (gitignored payload)
cisbp_database/         vendored CIS-BP MEME files (no external download needed)
docs/environment.md     full env-management guide (Tier 1/2/3: install, lockfile, conda-pack)
examples/demo.taiji.input.xlsx   sample of what build-taiji-input produces
templates/              SBATCH and Jinja templates for the agent's submit step
environment.yml         conda env spec (Python + R + Seurat/Signac + MACS3 + xlsx deps)
conda-lock.linux-64.yml lockfile for byte-reproducible installs (generate via `conda-lock`)
pyproject.toml          declares the `taiji-agent` package + `taiji-agent` CLI entry point
```

## Skills currently built

Eight production skills, all chainable end-to-end. The `taiji/` directory is reserved for a future umbrella skill but is empty today. For per-skill details — full input flags, output schemas, edge cases — read `skills/<name>/SKILL.md` and `skills/<name>/references/`.

### `sc-qc` — cell-level QC filtering

Applies per-cell quality filters to a Seurat object (RNA-only, ATAC-only, multiome same-cell, or co-embedded separate-cell) before it enters `coembed-construct` or `pseudobulk-construct`. This is the step the repo previously left entirely user-side.

Auto-detects modality from the object: `rna` (RNA assay only), `atac` (ATAC/peaks assay only), `multiome` (both assays, same cells — keeps intersection passing both filter sets), `coembed` (both assays, `meta.data$assay ∈ {RNA, ATAC}` — splits by origin, filters each separately, recombines).

**RNA filters:** `nFeature_RNA` strict min/max (> 200 & < 5000), `percent.mt` strict max (< 10%). `percent.mt` is auto-computed if absent. `nCount_RNA` is intentionally not filtered — matches the lab reference pipeline (`coembed_preprocess.R`). See `skills/CLAUDE.md` for the canonical threshold table.

**ATAC filters:** `nCount_<assay>` min, `nFeature_<assay>` min, `TSS.enrichment` min, `nucleosome_signal` max, `blacklist_ratio` max. Metrics absent from `meta.data` are **skipped with a warning** — they must be pre-computed upstream (Signac's `TSSEnrichment()`, `NucleosomeSignal()`, `FractionCountsInRegion()`).

Supports single-file mode (`--input`) and dual-file mode (`--rna-input` / `--atac-input`) for filtering separate RNA + ATAC objects before co-embedding.

Outputs: `filtered.rds` (or `rna_filtered.rds` + `atac_filtered.rds`), `qc_summary.json` (cells before/after per filter, marginal removal counts), `qc/violin_*_{before,after}.png`.

Explicit limitations: no doublet detection (run Scrublet/DoubletFinder upstream); fixed thresholds only (no MAD-based adaptive cutoffs).

### `detect-dataset-type` — pre-flight classifier

Classifies a directory as `bulk` / `single-cell` / `mixed` / `unknown` from file extensions plus 10x cellranger filename patterns. For single-cell, sub-classifies into `multiome` (RNA+ATAC same cells), `separate-assay` (different cells, needs Signac `integrate_atac`), or `sc-undetermined`. No file opens — fast and deterministic.

Recognized signatures: bulk via `.tsv` / `.narrowPeak` / `.bedpe`; single-cell via `.h5ad` / `.rds` plus 10x filename overrides (`fragments.tsv`, `barcodes.tsv`, `matrix.mtx`, `filtered_*_bc_matrix.h5`, etc.) which take priority over bulk extensions.

Used by `build-taiji-input` (refuses SC inputs) and `pseudobulk-construct` (refuses bulk inputs; uses `sc_modality` to pick clustering signal).

### `build-taiji-input` — bulk Taiji xlsx builder

Constructs the input xlsx (`Active` + `active_metadata` sheets) the bulk Taiji binary consumes. Walks `--data-dir` for per-sample RNA `.tsv`, ATAC `.narrowPeak`, optional HiC `.bedpe`; joins against a samples manifest; emits a properly-shaped xlsx.

Pre-flight gate soft-imports `detect-dataset-type` and refuses SC inputs with modality-specific guidance (multiome → sc-Taiji; separate-assay → coembed-construct; sc-undetermined → ask user to disambiguate). Degrades to warning if the sibling skill isn't installed.

**Three build-taiji-input bugs fixed (2026-04-30 small-test run):**
1. **`_peaks.rds` in `atac/` triggered the mixed-dataset detector** — `call_peaks.R` previously saved GRanges RDS files next to `.narrowPeak` files in `atac/`; the recursive `detect-dataset-type` scan saw `.rds` (SC) + `.narrowPeak`/`.tsv` (bulk) and refused. Fixed by: (a) saving `.rds` files to `atac/rds/` in `call_peaks.R`, and (b) switching `build-taiji-input`'s detect scan from `recursive=True` to `recursive=False` (top-level only avoids the `atac/rds/` subdir; the manifest `manifest.tsv` at the top level is still found as a `.tsv`/bulk file).
2. **Manifest TSV not parsed correctly** — `load_samples()` used `csv.DictReader` with default comma delimiter; `manifest.tsv` is tab-separated. Fixed by inferring delimiter from file extension (`.tsv`/`.tab` → `\t`, otherwise `,`).
3. **`--samples` path resolution** — when `--samples` is the pseudobulk `manifest.tsv`, pass `--rna-pattern "rna/{group}.tsv"` and `--atac-pattern "atac/{group}.narrowPeak"` so `build-taiji-input` resolves paths relative to `--data-dir`.

### `coembed-construct` — separate-assay scRNA + scATAC → shared latent space

When the user has two separate Seurat objects (one scRNA, one scATAC) for the same biological system, runs the [Stuart 2019 / Signac integrate_atac vignette](https://stuartlab.org/signac/articles/integrate_atac) end-to-end and emits a single coembed `.rds` that drops directly into `pseudobulk-construct` with `--use-existing-clusters`.

Pipeline: per-modality preprocessing (RNA NormalizeData→PCA→UMAP; ATAC TFIDF→SVD→UMAP) → ATAC `GeneActivity` → `FindTransferAnchors(reduction="cca")` → `TransferData` of RNA matrix to ATAC cells → `merge()` with `meta.data$assay ∈ {RNA, ATAC}` origin tag → joint PCA/UMAP → cluster.

What it deliberately does NOT do: label transfer onto cells. Cluster IDs are de novo; named cell-type labels are a user precondition (judgment-required upstream QC stays user-side).

Outputs: `coembed.rds` (drop-in for pseudobulk-construct), `coembed_summary.json` (anchor count, feature-drop audit, peak-RAM estimate, cluster table), `qc/{umap.png, umap_coords.csv}`.

Hardened against several silent-failure classes — see `skills/coembed-construct/SKILL.md` and references/ for: metadata-value casing mismatches across modalities, Seurat v5 inputs in v4 envs, stale ChromatinAssay fragment paths, macOS R `R_MAX_VSIZE` cap, macOS jetsam OOM under swap pressure, FindTransferAnchors feature-drop audit, Louvain modularity saturation, and the GenomeInfoDb / biovizBase Bioconductor deps that aren't always re-exported under suppressPackageStartupMessages.

**HPC note — `seqlevelsStyle` network call:** `seqlevelsStyle(annotations) <- "UCSC"` triggers a UCSC chromInfo download; on HPC nodes without outbound internet this fails. `coembed.R` wraps this in a `tryCatch` with a manual `chr`-prefix fallback (MT → chrM), so HPC runs proceed without internet. If you see `cannot open URL 'https://hgdownload.soe.ucsc.edu/...'` in output, this fallback is working correctly.

**`--output` path must include `.rds`:** `coembed.py` now auto-appends `.rds` if the supplied path has no extension, so `--output runs/name/coembed` becomes `runs/name/coembed.rds`. Pass the full path with extension to avoid the NOTE message.

### `pseudobulk-construct` — single-cell → bulk-Taiji bridge

Converts a single-cell object (`.rds` / `.h5ad`) plus a fragments file into pseudobulk RNA GeneQuant TSVs and per-cluster narrowPeak files in the exact layout `build-taiji-input` consumes. Closes the SC → bulk Taiji loop.

Pipeline (orchestrated by `pseudobulk.py`): detect-dataset-type gate → dependency check → `load_and_cluster.R` (Seurat/Signac WNN clustering or RNA-only PCA / ATAC-only LSI; scale-aware resolution binary search targeting mean ~200 cells/cluster, drops <20-cell clusters) → `aggregate_rna.R` (sums raw counts per (cluster × metadata) group → 2-column GeneQuant TSV) → `call_peaks.R` (per-cluster `Signac::CallPeaks`; barcode reconciliation against fragments.tsv.gz with ≥95% overlap required, fail-loud otherwise; emits `<group>.narrowPeak`) → `manifest.tsv` for direct hand-off to build-taiji-input.

**`call_peaks.R` — per-group fragment extraction:** `call_peaks.R` uses a `zcat|awk` pre-filter to write a cell-specific temp BED per group before calling MACS3 directly (bypassing `Signac::CallPeaks`). This produces per-group peaks and is ≈10–100× faster on large fragment files than using `CallPeaks` which passes the full file path to MACS3. `Signac::CallPeaks` ignores the `cells=` argument when writing the MACS3 input — it passes the raw fragment file, so WITHOUT this pre-filter all groups would share the same peak set.

Stratification: `--metadata-cols` cross-products within each cluster (e.g. `genotype,tissue` with values {WT, KO} × {site_a, site_b} → four per-cluster groups). `--cohort-col` picks which metadata column becomes the manifest cohort label.

Output layout under `<output-dir>/`: `clusters.csv`, `resolution_trace.json`, `groups_plan.json`, `rna/<name>.tsv`, `atac/<name>.narrowPeak`, `atac/rds/<name>_peaks.rds`, `qc/umap.png`, `qc/umap_coords.csv`, `manifest.tsv`. Note: GRanges peak objects go in `atac/rds/` (not directly in `atac/`) so that `build-taiji-input`'s `detect-dataset-type` scan sees a clean bulk directory.

### `workflow-log` — per-run audit trail

Constructs and updates a per-run log (`.md` + `.jsonl`) under `<work_dir>/log/` capturing every stage of a Taiji-agent workflow. Sibling skills auto-attach via a small `active_run` pointer file — when a run is active, every call to `detect`, `pseudobulk`, `build-taiji-input`, `fetch-references`, or per-sample `taiji-runner` writes a stage entry; when no run is active, the skills behave exactly as before (zero-cost soft import).

Outputs per run in `<work_dir>/log/`: `<run_id>.md` (narrative), `<run_id>.jsonl` (machine-parseable), `index.jsonl` (one line per run for cross-run grep), `active_run` (pointer; deleted on finalize). `<run_id>` = `YYYYMMDD_HHMMSS_<4hex>` (sortable, collision-resistant).

```python
from log import TaijiLog
log = TaijiLog.start(work_dir=".", genome="hg38",
                     reference_paths={"fasta": "...", "gtf": "...", "meme": "..."})
log = TaijiLog.attach_active(work_dir=".")     # returns None if no active run
if log: log.append_detect(detect_result)
log.finalize(status="success", duration_s=...)
```

Soft-attach contract: sibling skills always check for `None` before calling `append_*`; logging failures inside the sibling skill are swallowed so a log bug never breaks a workflow.

### `fetch-references` — per-genome reference-data stager

Idempotently stages the reference data Taiji needs into `<output_dir>/dependencies_data/`. FASTA + GTF are downloaded from GENCODE via the EMBL-EBI mirror; the MEME motif file is **CIS-BP** (vendored at `cisbp_database/`, no network). Reads a single `reference_manifest.yml` declaring URLs, target paths, gunzip flags, and approximate uncompressed sizes per genome. **HOCOMOCO is not supported** — locking the project to a single motif source keeps PageRank values comparable across runs.

Genomes available out of the box: `hg38` (GENCODE v45), `hg19` (GENCODE v19), `mm10` (GENCODE M25). Add new ones by editing `scripts/reference_manifest.yml`.

Idempotency: files at the expected path with size within ±5% of the manifest's `approx_size_mb` are skipped. Downloads are atomic (`.part` + rename). `--force` re-downloads; `--check` reports without downloading; `--dry-run` resolves URLs without touching the network. The skill never deletes files it didn't write itself.

`--update-genomes-yml` rewrites `skills/build-taiji-input/assets/genomes.yml` in place to point at the resolved fasta/gtf paths. `samtools faidx` is auto-built if samtools is on PATH.

### `taiji-runner` — per-sample Taiji 1.3.0 orchestrator

Bridges between `build-taiji-input`'s xlsx (the human manifest) and Taiji's actual input format. `taiji run --config <yml>` is single-sample-per-call and reads YAML/TSV, not xlsx — so a multi-sample run requires splitting the xlsx into N TSVs and N configs, then invoking Taiji N times. This skill does that.

Mirrors the existing UCSD wrapper at `Taiji/taiji_wrapper-uniqueGen.py` — same directory layout (`Input/<sample>_input/`, `Output/Partial/<sample>_output/`), same three-placeholder convention (`[insert_input_filepath_here]`, `[insert_output_directory_here]`, `[insert_genome_filepath_here]`), same `taiji_config_files.txt` index. Drop-in compatible.

Output validation against silent failures: Taiji 1.3.0 sometimes exits 0 even when its workflow errored. The runner cross-checks: a sample is `ok` only if `exit_code == 0` AND `GeneRanks.tsv` exists in its output dir. Otherwise the run is flagged as a silent workflow failure. **Critical naming gotcha**: Taiji 1.3.0 emits `GeneRanks.tsv`, NOT `*_pagerank.tsv`. The binary's symbol table contains stage names like `Output_Ranks` / `Output_Ranks_SC`, but those are workflow STAGES — the actual output file is `GeneRanks.tsv`.

Verified Taiji 1.3.0 output schema (`Output/Partial/<sample>_output/`): `GeneRanks.tsv` (per-TF PageRank — ~1100 rows for cisbp_human_2.meme on hg38), `GeneRanks_PValues.tsv`, `sciflow.db` (workflow state DB), `Network/<sample>/{edges_combined.csv, edges_binding.csv, nodes.csv}` (Neo4j-style CSVs), `ATACSeq/{openChromatin.bed.gz, TFBS/}`, `RNASeq/{expression_profile.tsv, RNA-<sample>_rep1_gene_quant.tsv}`, `GENOME/genome.index` (~3 GB; built once per genome, reused across samples).

Composition with `bin/run-taiji.sh`: the shell wrapper picks the binary by OS, regenerates the xlsx, runs preflight, starts/finalizes the workflow log, and delegates per-sample to taiji-runner. Users who want only the per-sample machinery can call `python skills/taiji-runner/scripts/run_taiji.py` directly.

Executor auto-detection: `--executor auto` (default) picks `slurm` when `sbatch` is on PATH, else `local`. SLURM mode submits an array job (`sbatch -a 1-N%MAX`, default MAX=10), polls `squeue` until the array completes, then validates outputs the same way as local mode. `--no-wait` returns the job ID immediately. Local mode stays sequential by default for memory safety; SLURM mode defaults to `--parallel 10`.

## Run-directory layout

Per-run convention (also documented in skills/build-taiji-input/SKILL.md and skills/taiji-runner/SKILL.md): everything a human authors for a run lives under `runs/<name>/code/` — the `samples.csv` manifest, the per-sample `taiji_config.template.yml`, plus any ad-hoc preprocessing or post-run analysis scripts. Anything outside `code/` is generated by the skills (Input/, Output/, log/, taiji_input.xlsx, taiji_config_files.txt) and is gitignored. `bin/run-taiji.sh` resolves `samples.csv` and `taiji_config.template.yml` from `<run_dir>/code/` first; the flat layout (`<run_dir>/samples.csv`) is still accepted for back-compat.

## Skill composition (end-to-end pipeline)

```
    fetch-references --genome <build> --output dependencies_data/  (one-time per site)
                              │
                              ▼
    workflow-log.start(work_dir, genome, reference_paths) ──┐
                              │                              │
                       detect-dataset-type ─┐                │
                              │             ├── auto-attach to active run
            ┌─────────────────┴────────┐    │  via <wd>/log/active_run
            ▼                          ▼    │  pointer; emit stage entries
      classification=bulk      classification=single-cell
            │                          │
            │           ┌──────────────┴──────────────┐
            │           ▼                             ▼
            │    sc_modality=multiome     sc_modality=separate-assay
            │           │                             │
            │         sc-qc                  sc-qc (--rna-input / --atac-input)
            │      (--input obj.rds)         filter each modality separately
            │           │                             │
            │           │                  coembed-construct ──┐
            │           │                  (Stuart/Signac      │
            │           │                   integrate_atac)    │
            │           │                             │        │
            │           │                       coembed.rds    │
            │           │                             │        │
            │           ▼                             ▼        │
            │    pseudobulk-construct  ◄───  pseudobulk-construct
            │    (--clustering-signal wnn)   (--use-existing-clusters)
            │           │                             │        │
            │           └──────────────┬──────────────┘        │
            │                          │                       │
            │                       manifest.tsv               │
            │                          │                       │
            └────────►   build-taiji-input  ◄──────────────────┤
                              │                                │
                              ▼                                │
                  taiji_input.xlsx (Active + active_metadata)  │
                              │                                │
                              ▼                                │
              taiji-runner: split xlsx → per-sample TSVs/configs ─┤
                              │                                   │
                              ▼                                   │
              for each sample: bulk Taiji binary                  │
                              │                                   │
                              ▼                                   │
              workflow-log.finalize(status) ──────────────────────┘
```

Each step gates on the previous one and produces inputs in the next step's expected format with **no glue code in between**. The `manifest.tsv` from `pseudobulk-construct` is a drop-in for `build-taiji-input --samples`.

## External dependencies (not vendored)

All inputs to Taiji are pre-aligned/pre-quantified — no BWA, STAR, RSEM, samtools, Picard, or Singularity needed. Per-skill dep manifests live at `skills/*/dependencies.yml`; `bash bin/doctor.sh` verifies them.

**Always required:**
- Python ≥ 3.10
- Pinned Python deps in `environment.yml`: click, pydantic, pyyaml, jinja2, rich, pandas, openpyxl
- Rscript ≥ 4.2 (only for `pseudobulk-construct` and `coembed-construct`)
- R packages (bioconda): Seurat ≥ 5.0, Signac ≥ 1.12, Matrix, GenomicRanges, GenomeInfoDb, biovizBase, dplyr, optparse, jsonlite
- SeuratDisk (for `.h5ad`) is **not** auto-installed — unmaintained upstream, brittle to build on cluster envs. `.h5ad` ingestion requires either a manual `remotes::install_github("mojaveazure/seurat-disk")` or upstream conversion to `.rds`.
- MACS3 ≥ 3.0 (or MACS2 — the SC skills auto-detect whichever is on PATH)
- gunzip + awk (standard Unix)

**Taiji binary** — installed by `bin/install-taiji.sh --system <centos|ubuntu|macos>`. Singularity is explicitly out of scope (no official .sif published).

## Environment management

Composable profiles let users install only what their dataset needs. Bulk-only users save ~3.5 GB and ~20 min by skipping the R + Seurat/Signac stack.

```bash
bash bin/install.sh --system <macos|centos|ubuntu>       # base only (default)
bash bin/install.sh --system macos --profile sc          # base + sc (additive)
bash bin/auto-install.sh --data-dir <path> --system macos        # detect-driven
micromamba activate taiji-agent
bash bin/doctor.sh --profile base                        # filter by profile
```

| Profile | Contents | Enables | Disk | Time |
|---------|----------|---------|------|------|
| `base` (default) | Python + click/pydantic/pyyaml + pandas + openpyxl + macs3 + local pkg | detect-dataset-type, build-taiji-input, fetch-references, taiji-runner, workflow-log | ~500 MB | ~5 min |
| `sc` (additive) | r-base + r-seurat + r-signac + Bioconductor (GenomicRanges + GenomeInfoDb + biovizBase + EnsDb) | sc-qc, pseudobulk-construct, coembed-construct | +3-4 GB | +15-30 min |
| `dev` (orthogonal) | pytest + pytest-cov + ruff + mypy + ipython | author tooling | +500 MB | +3 min |

`sc` is **additive** — a user who installed `base` and later needs SC can run `bash bin/install.sh --profile sc` and only the R packages get added. Each skill declares its profile in `skills/<name>/dependencies.yml`; `bin/doctor.sh --profile <name>` filters the verification table.

Reproducibility tiers (full details in `docs/environment.md`): Tier 1 micromamba + environment.<profile>.yml; Tier 2 conda-lock for byte-identical envs across machines; Tier 3 conda-pack tarball for a portable, network-free install.

## Design conventions across all skills

These are load-bearing — preserve them when adding new skills or extending existing ones.

- **CLI is a thin wrapper over a pure library function.** Every skill exposes a callable Python (or R) function returning a dataclass; the CLI is just argparse + `print(result.summary())`.
- **No file content sniffing.** All classification and routing happens on extensions and filenames. The cost of opening files (h5py, R runtime) is not worth the marginal robustness.
- **Strict-on-mixed by default.** Silent SC contamination of a bulk pipeline produces garbage at Taiji-invocation time. Loud-fail-early > debug-at-analysis-time.
- **Soft-import siblings.** Cross-skill dependencies use `sys.path.insert(0, sibling/scripts)` and a `try: from <sibling> import ...` so a standalone install of one skill still works without the other.
- **Reference docs in `references/`, not in `SKILL.md`.** SKILL.md stays short; deep edge-case reasoning lives in `references/<topic>.md`.
- **Every recognized format/pattern ships with a fixture and an assertion.** `<skill>-workspace/fixtures/` mirrors `evals/evals.json`. When a verification pass finds a bug, add the failing case as a fixture.
- **Always cite limitations.** Explicit "this heuristic can fail when X" notes over silent confidence.
- **Output formats designed for downstream chaining.** When skill A's output feeds skill B, A emits exactly what B expects with no transformation step in between.

## Single-cell preprocessing conventions (project-wide)

These apply to every `.rds` / `.h5ad` that flows into `coembed-construct` or `pseudobulk-construct`. The skills detect violations at startup and warn, but harmonizing upstream is strictly better than relying on the runtime warning — it eliminates the failure class entirely instead of just surfacing it.

- **Lowercase every metadata-value string before saving.** Casing mismatches across modalities (e.g. `Tissue_A` on RNA vs `tissue_a` on ATAC) survive `merge()` as TWO distinct values and silently halve effective sample size in cross-product stratification — the skills emit `WARN: CASE-ONLY mismatch`, but the user-side fix is what removes the bug. Idiomatic R applied to BOTH per-modality objects before they enter the pipeline:
  ```r
  for (col in c("tissue", "genotype", "sample_id", "donor", "batch")) {
    if (col %in% colnames(obj@meta.data)) {
      obj[[col]][[1]] <- tolower(as.character(obj[[col]][[1]]))
    }
  }
  saveRDS(obj, "obj_clean.rds")
  ```
  This also collapses `WT`/`wt`, `Day7`/`day7`, etc. Free-text columns not passed through `--metadata-cols` are unaffected.
- **Use the same metadata column NAME for the same biological concept across modalities.** The skill validators only see what you ask them to compare — `sample_id` on RNA + `library_id` on ATAC for the same biological samples is invisible.
- **Rebuild stale ChromatinAssay fragment paths before shipping `.rds` across machines.** `Fragments(atac)` embeds an absolute path at object-build time. The `coembed-construct` skill accepts `--fragments` as a runtime override; the cleaner upstream fix is to rebuild the Fragment handle on the destination machine before saving.
- **Match the Seurat version of the input `.rds` to the env that consumes it.** Saving under v5 then loading under v4 silently produces an empty `Assays()` and crashes the pipeline ~10 min in. The skill refuses up-front under v4-with-v5-input; the user-side fix is to upgrade the env or downgrade the input via `obj[['RNA']] <- as(obj[['RNA']], 'Assay')`.

## Memory system

Persistent memory at `~/.claude/projects/.../memory/` (path is environment-specific; check `MEMORY.md` in your session). Contains user role and preferences, project decisions, skill-specific design notes, and feedback memories. Read these before making non-trivial decisions about extending the project.

## Sandbox model (workspace-bound execution)

Every command the agent invokes runs through `bin/sandbox-run.sh`, a layered wrapper that ensures all writes stay inside the workspace folder. Three tiers, applied in order: (1) **soft** (always): cwd must be inside `REPO_ROOT`, `TMPDIR` pinned to `<REPO_ROOT>/tmp/sandbox-<pid>`; (2) **OS-level** (auto): `bwrap` on Linux or `sandbox-exec` on macOS, restricting fs writes to `REPO_ROOT` + workspace TMPDIR and blocking network unless `--allow-net`; (3) **strict** (opt-in): `--strict-sandbox` or `TAIJI_STRICT_SANDBOX=1` refuses to run unless an OS-level sandbox tool is present.

`bin/run-taiji.sh` already wraps every Taiji invocation with `sandbox-run.sh`. The Taiji `tmp_dir` is pinned to `${REPO_ROOT}/tmp/<run_name>` so even Taiji's scratch I/O stays in the workspace. `--no-sandbox` exists for debugging; `--strict-sandbox` for shared infrastructure where soft enforcement isn't enough. Threat model: agent mistakes and buggy scripts, NOT malicious-code escape — for the latter, layer Docker/Apptainer on top.

## Architectural constraints (load-bearing)

- **The Taiji binary cannot be executed from a sandbox host with a different ISA.** When the host is aarch64 (e.g. an Apple-Silicon-based agent sandbox) and binaries are x86_64 ELF/Mach-O with no qemu-user / box64 / fex-emu emulator installed, the kernel rejects with `Exec format error`. Practical consequence: every "run Taiji" step is either Mac-side or SLURM-side. The agent's role from such a sandbox is strictly **prepare inputs + validate outputs**.
- **Workspace folder is shared between the sandbox and the host machine.** Files dropped into `<workspace>/data/`, `<workspace>/dependencies_data/`, `<workspace>/runs/<name>/Output/` from either side become immediately visible to the other. Use this for handoffs: host runs Taiji → outputs land in workspace → sandbox can validate / summarize / log without re-running.
- **R may be available in the sandbox without the SC stack.** Base R `readRDS` lets you inspect Seurat object slots / class / metadata columns / barcode formats WITHOUT installing the full Seurat/Signac/Bioconductor stack — useful for pre-flight diagnosis (catching multi-sample barcode suffixes, pre-existing reductions, unexpected assay layouts) before submitting a long SLURM run. Full coembed/pseudobulk pipelines belong on Mac/SLURM.
- **Seurat v5 API breaking changes in R scripts.** `GetAssayData(slot=)` is defunct in SeuratObject ≥5.0.0 — use `layer=` (wrap in `tryCatch` to stay v4-compatible). `seqlevelsStyle(x) <- "UCSC"` in GenomeInfoDb triggers a network fetch from UCSC; HPC nodes without internet need a manual fallback (`seqlevels(x) <- ifelse(grepl("^chr", lvls), lvls, ifelse(lvls=="MT","chrM",paste0("chr",lvls)))`). R does NOT support implicit string concatenation — multi-line help text in `make_option()` must be a single string or use `paste0()`.
- **`bgzip`/`tabix` are not in the taiji-agent conda env.** When pre-filtering fragments files, use `/stg3/data1/eunice/bin/miniconda3/envs/bcftools/bin/bgzip` and `../tabix`. Signac's `CreateFragmentObject` requires the `.tbi` index to exist but does not need bgzip on PATH.

## Not yet built (scoped but pending)

- The actual Taiji CLI plugin wrapping these skills (the `skills/taiji/` directory is reserved for this).
- A SLURM submission skill that takes a `taiji_input.xlsx` + a config and emits an sbatch script.
- A cron watchdog that monitors upstream output directories (e.g. cellranger output dirs) and triggers the pipeline when new datasets land.
- For `pseudobulk-construct`: real fixture `.rds` / `.h5ad` files, an eval loop matching `build-taiji-input`'s and `detect-dataset-type`'s, and first-run validation on a SLURM env (R syntax check + end-to-end on a public PBMC multiome).

## Where to find what

| Question                                  | Where                                                                 |
|-------------------------------------------|-----------------------------------------------------------------------|
| "What does skill X do?"                   | `skills/X/SKILL.md`                                                   |
| "Why does skill X make decision Y?"       | `skills/X/references/`                                                |
| "What design decisions were locked in?"   | persistent memory (see Memory system section)                         |
| "What fixtures cover what cases?"         | `skills/X-workspace/fixtures/` + `skills/X/evals/evals.json`          |
| "How is the Python package wired up?"     | `pyproject.toml` + `src/taiji_agent/`                                 |
| "What conda env is expected?"             | `environment.yml`                                                     |
| "What does a Taiji input xlsx look like?" | `examples/demo.taiji.input.xlsx`                                      |
