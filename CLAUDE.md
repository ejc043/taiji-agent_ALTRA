# Taiji agent — repository guide for Claude

## What this repo is

A Python package + plugin-style skills bundle that wraps the [Taiji pipeline](https://github.com/Taiji-pipeline/Taiji) (Wei Wang lab, UCSD) for routine use on a SLURM HPC. The agent does not re-implement Taiji — it handles preflight validation, input construction, single-cell preprocessing, and submission so the user can hand it a directory of data and get a Taiji-ready run with minimal manual setup.

User: Eunice Chen (UCSD bioinformatician, ejc043@ucsd.edu). Practitioner-level — appreciates context, summaries, and explicit limitations. Not a beginner coder.

## Top-level layout

```
src/taiji_agent/        Python package (CLI scaffold; cli, config, detect, samplesheet, validate)
skills/                 Plugin-style skills (the bulk of the project's logic lives here)
  build-taiji-input/        produces the bulk-Taiji input xlsx
  coembed-construct/        co-embeds separate-assay scRNA + scATAC into shared latent space (Stuart/Signac integrate_atac vignette)
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
  postinstall.R             remotes::install_github() for MuDataSeurat (.h5mu input)
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

Seven production skills, all chainable end-to-end. The `taiji/` directory is reserved for a future umbrella skill but is empty today.

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

### 3a. `coembed-construct` — separate-assay scRNA + scATAC → shared latent space

**Purpose:** When the user has two separate Seurat objects (one scRNA, one scATAC) for the same biological system, run the [Stuart 2019 / Signac integrate_atac vignette](https://stuartlab.org/signac/articles/integrate_atac) end-to-end and emit a single coembed `.rds` that drops directly into `pseudobulk-construct` with `--use-existing-clusters`.

**Pipeline (`coembed.R`):**
1. RNA: NormalizeData → FindVariableFeatures → ScaleData → RunPCA → RunUMAP.
2. ATAC: set EnsDb annotation (hg38/hg19 → `EnsDb.Hsapiens.v86`; mm10/mm39 → `EnsDb.Mmusculus.v79`) → RunTFIDF → FindTopFeatures → RunSVD → RunUMAP.
3. `GeneActivity` → ACTIVITY assay on ATAC, normalize + scale.
4. `FindTransferAnchors(reference = rna, query = atac, query.assay = "ACTIVITY", reduction = "cca")`.
5. `TransferData` of the RNA expression matrix → imputed RNA assay on ATAC cells (using the `lsi` weight reduction).
6. `merge(rna, atac)` → coembed object; `meta.data$assay` ∈ {"RNA", "ATAC"} tags origin.
7. Joint ScaleData → RunPCA → RunUMAP on RNA features (real for RNA-origin cells, imputed for ATAC-origin).
8. FindNeighbors + FindClusters with the same scale-aware resolution binary search `pseudobulk-construct` uses (target ~200 cells/cluster, drops <20-cell clusters). Override with `--resolution <r>` or skip with `--no-cluster`.

**What it deliberately does NOT do:** label transfer onto cells (`TransferData(refdata = seurat_annotations)`). Cluster IDs come from de novo clustering on the shared space; named cell-type labels are user-precondition (consistent with the "judgment-required upstream QC stays user-side" memory).

**Outputs:**
- `<output_dir>/coembed.rds` — drop-in for `pseudobulk-construct --input <path> --use-existing-clusters --fragments ...`
- `<output_dir>/coembed_summary.json` — anchor count, gene-list size, chosen resolution, cluster table by `(cluster × assay)`.
- `<output_dir>/qc/{umap.png, umap_coords.csv}` — UMAP panels colored by clusters / assay / metadata cols.

**Status:** R script syntax-checked. R 4.2.0 is available in the sandbox at `/usr/local/bin/Rscript`, but Seurat / Signac / EnsDb are NOT pre-installed (slow Bioconductor compile from source — ~30 min). For real coembed runs, install on the target machine via `bin/install.sh --profile sc` and run there. Eunice's SLURM env has the full SC stack.

**Real-data gotchas surfaced from a 60k+21k cell run** (D7 mouse, 4-sample-aggregated barcodes with `-1`/`-2`/`-3`/`-4` suffixes):
- Inputs often arrive with their own pre-computed `pca`/`umap` (RNA) and `lsi`/`umap` (ATAC) reductions PLUS per-modality `seurat_clusters`. Recomputing from scratch is the safe default but wastes ~10 min — pass `--reuse-rna-reductions` and/or `--reuse-atac-reductions` when the inputs are clean. The skill always preserves the input cluster labels as `rna_input_clusters` / `atac_input_clusters` (non-clobbering) so the merged object lets you compare per-modality vs joint clustering downstream.
- Multi-sample-aggregated barcodes (10x Multi / cellranger-aggr) have suffixes `-1`/`-2`/.../`-N` matching across the RNA object, ATAC object, and fragments file. Identity-match works as-is; the skill's suffix-stripping reconciliation logic doesn't kick in (correct).
- ATAC objects sometimes arrive with a pre-existing `RNA` assay (from a prior integrate_atac, or sliced multiome). The skill OVERWRITES it with fresh TransferData-imputed values and logs a NOTE — if you want to preserve the original, copy it under a different assay name first.
- **Metadata-value casing mismatches across inputs are the highest-priority silent bug.** Real example from this dataset: RNA has `tissue ∈ {spleen, siIEL}`, ATAC has `tissue ∈ {Spleen, siIEL}`. After merge, `tissue` has THREE distinct values; cross-product stratification in `pseudobulk-construct` then creates separate groups for `spleen` and `Spleen`, halving the effective sample size with neither group flagged as wrong. The skill now validates value sets across both inputs at startup with `validate_metadata_consistency()` and logs `WARN: CASE-ONLY mismatch` (or `WARN: VALUE mismatch`) plus a copy-paste fix hint. Use `--strict-metadata` for production runs to convert this into a hard refusal.
- **Seurat v5 (Assay5) inputs are supported but look empty in raw probes.** `attributes(obj)$assays$RNA$counts` is empty on v5 — the data lives in `attributes(.)$layers$counts`. `GetAssayData()` dispatches across v3/v5 in Seurat ≥ 5.0, so the skill works unchanged in a proper v5 env. Mixed-version inputs (v5 RNA + v3 ATAC) merge fine in Seurat ≥ 5.0; if `merge()` errors, run `JoinLayers(rna)` upstream.
- **The hardest silent v5 failure: v5 inputs in a v4 env.** Verified on the D7 dataset against Seurat 4.3.0: `readRDS` succeeds, `class(obj)` is "Seurat", but `Assays(obj)` returns empty and `GetAssayData(obj, slot="counts")` errors with `'RNA' is not an assay` ~10 minutes into the pipeline. The skill now checks `packageVersion("Seurat")` and the per-assay class at startup and refuses with a clear `Either upgrade to Seurat >= 5 or downgrade the input via obj[['RNA']] <- as(obj[['RNA']], 'Assay')` message — saving the user 30 min of wasted compute. Common triggers: lab member saved a .rds in a Seurat v5 env then a different lab member tries to consume it in their older anaconda env that pins to Seurat 4.x. The bioconda recipe in `environment.sc.yml` pins to `r-seurat>=5.0` to prevent this in the canonical env.
- Joint PCA on >50k merged cells needs 10–20 GB peak RAM. SBATCH `--mem=128G` is recommended for >80k merged cells. See `references/coembed_strategy.md` for the resource-estimation reasoning.
- **Stale fragments path inside ChromatinAssay (HPC → laptop transfer wart).** The D7 ATAC `.rds` was built on `/home/whw006/...`; on Eunice's Mac that path doesn't exist, but `readRDS` succeeds silently. `GeneActivity()` then fails ~10 min in with an opaque tabix error. The skill now `file.exists()`-checks every Fragment object's `path` slot at startup and accepts `--fragments` to rebuild a single Fragment handle in place (reusing the existing barcode-cell mapping). If `--fragments` is not provided and a path is missing, the skill refuses loudly with a copy-paste fix hint. Common trigger: `.rds` shipped from a collaborator's HPC; canonical fix is to call this skill with `--fragments <path>` rather than mutate the input `.rds`.
- **`GetAssayData(slot = ...)` is defunct in SeuratObject ≥ 5.0** (not just deprecated — hard error via `lifecycle::deprecate_stop`). Skills must use `GetAssayData(obj, layer = "counts" | "data" | "scale.data")`. The slot/layer migration was completed in coembed-construct after the D7 dataset triggered the lifecycle stop on the very first GetAssayData call.
- **Hidden Bioconductor deps for Signac annotation extraction.** `Signac::GetGRangesFromEnsDb()` calls into `biovizBase::crunch()`, and `seqlevelsStyle(annotations) <- "UCSC"` needs `GenomeInfoDb` (its assignment generic is not always re-exported under `suppressPackageStartupMessages`). Both are now in `skills/coembed-construct/dependencies.yml` and the R script `library()`s `GenomeInfoDb` explicitly. Without these the pipeline dies at the EnsDb load step right before GeneActivity.
- **macOS R has a default `R_MAX_VSIZE` cap (~16-18 GB) below what `FindTransferAnchors` needs on a 60k+21k pair.** Linux R has no such cap. R errors with `vector memory limit of 18.0 Gb reached, see mem.maxVSize()` partway through CCA. The skill bumps the cap to 128 GB at script startup on Darwin only (override via `TAIJI_R_MAX_VSIZE_GB`). On 18 GB-physical Macs the spike lands in macOS swap and the run *can* complete — but only if the swap budget isn't already eaten by other apps (see jetsam-OOM gotcha below). SLURM/Linux nodes are unaffected.
- **macOS jetsam OOM-kills R mid-CCA when swap budget is exhausted.** This is distinct from the VSIZE cap above and worse: bumping VSIZE prevents R from self-throttling, but the kernel can still kill R via SIGKILL/SIGTERM when **physical RAM + swap budget** are both saturated. Verified on the D7 dataset on Eunice's 18 GB M-series Mac on 2026-04-28: FindTransferAnchors peaked at ~10 GB RSS, swap filled to 23.7/24.6 GB (97%), the R process went into uninterruptible-sleep (state `U`) thrashing at ~3% CPU. Process exited with `-15` (SIGTERM) and **NO R error message** — the user just sees `coembed.R exited -15` after 5+ wasted minutes. The skill now (a) computes a pessimistic peak-RSS estimate `8 * n_features * (n_rna + n_atac) * 6 bytes` after both objects load, (b) sysctls `hw.memsize` + `vm.swapusage` for the actual macOS budget, (c) emits `INFO`/`WARN`/`HARD-WARN` based on `peak vs RAM + 0.5*free_swap`, and (d) writes the estimate to `coembed_summary.json` as `peak_ram_estimate_gb` for post-run audit. Pass `--abort-on-memory-risk` to refuse instead of warning. The Python wrapper recognizes signal-driven exits (-9 / -15 / 137 / 143) at the Rscript layer and prints a multi-line OOM diagnostic with concrete mitigations: run on SLURM (`sbatch --mem=128G ...`), `subset(rna, downsample=20000)` (the RNA reference dominates the CCA working set; ATAC stays full-size), lower `--n-pcs` from 30 to 20, free swap by closing GUI apps. In short: 18 GB Macs CAN run 60k+21k, but only when not memory-constrained by other workloads — for reliability, use SLURM.
- **High `FindTransferAnchors` feature drop rate signals an annotation mismatch.** `FindTransferAnchors(features = VariableFeatures(rna), ...)` silently drops every feature absent from EITHER assay. Because `Signac::GeneActivity` only computes coverage for genes with EnsDb annotations, baseline drop is ~10-20% (verified D7: 1639/2000 = 82% kept on mm10). Drop > 50% almost always means: wrong `--genome`, EnsDb species mismatch (mm10 EnsDb on human RNA), or gene-name convention mismatch (Ensembl IDs vs gene symbols on the RNA side). The skill now computes the intersection up-front, logs `features for anchoring: N requested -> M shared (X% kept, Y dropped)`, emits a `WARN` with the most likely root causes when drop > 50%, and records counts in `coembed_summary.json` as `n_features_requested` / `n_features_used` / `n_features_dropped`.
- **Louvain modularity saturation on large merged objects.** On 81k+ cells with `--target-cluster-size 200`, Louvain hits a modularity ceiling around 75-100 clusters and refuses to split further. The binary search now early-stops when `n_clusters` is unchanged across two iters at the upper bracket (`r=hi=5.0`), saving 5-7 min of identical FindClusters calls. Verified on the D7 dataset: 8 iters → 3 effective iters before saturation detection. If the user genuinely wants finer clusters, sub-cluster post hoc with `FindSubCluster()`.
- **Verified end-to-end on D7 (60k RNA + 21k ATAC, mm10).** Run took ~25-30 min on an 18 GB-physical M-series Mac (swap-bound during CCA), produced 67 088 anchors, 75 raw clusters / 74 retained at saturation res=5.0. Joint UMAP shows good RNA/ATAC mixing across most of the embedding space. The case-only `tissue` mismatch (`spleen` vs `Spleen`) created the predicted 3-tier split visible in `qc/umap.png` — surfacing this pre-run via `validate_metadata_consistency` is exactly what saved a future `pseudobulk-construct --metadata-cols tissue` from silently doubling the spleen sample count. A subsequent re-run on 2026-04-28 with `--reuse-rna-reductions --reuse-atac-reductions` skipped the per-modality preprocessing successfully but jetsam-killed during CCA when other apps were eating swap — leading to the new memory pre-flight + OOM-diagnostic gotchas above.
- **Verified again with the hardened skill on 2026-04-28.** Re-ran on a downsampled 15k+6k subset of D7 (to fit the same 18 GB Mac with no swap headroom) and confirmed the new pre-flight + audit logs fire correctly: peak-RAM estimate `1.9 GB at 21000 cells x 2000 features`, macOS budget `18 GB RAM + 1.5 GB free swap = 18.8 GB usable` (no warning since 1.9 < 18.8), features-for-anchoring audit `2000 requested -> 1639 shared (82% kept, 361 dropped)`, plus all the previously-handled gotchas (stale fragments rebuilt via `--fragments`, pre-existing RNA assay flagged, pre-existing seurat_clusters preserved as `*_input_clusters`, CASE-ONLY tissue mismatch warned with copy-paste fix). End-to-end run time: ~3.5 min on the downsample. `coembed_summary.json` now carries `n_features_requested` / `n_features_used` / `n_features_dropped` / `peak_ram_estimate_gb` for post-run audit. The 2×2 joint UMAP visibly contained 3 tissue colors (`siIEL` + `spleen` + `Spleen`), proving the validate_metadata_consistency warning identifies exactly the silent-bug case it claims to.

### 3. `pseudobulk-construct` — single-cell → bulk-Taiji bridge

**Purpose:** Convert a single-cell object (`.rds`, `.h5ad`, `.h5mu`) plus a fragments file into pseudobulk RNA GeneQuant TSVs and per-cluster narrowPeak files in the exact layout `build-taiji-input` consumes. Closes the SC → bulk Taiji loop.

**Inputs:**
- `--input` path to single-cell object.
- `--fragments` path to `fragments.tsv.gz` (required for ATAC peak calling).
- `--genome` tag.
- `--output-dir`.
- Optional: `--metadata-cols` (cross-product within each cluster — e.g. `genotype,tissue` → per-cluster `(WT,spleen)`, `(WT,siiel)`, `(KO,spleen)`, `(KO,siiel)` quads), `--cohort-col` (which metadata col's value becomes the manifest cohort label; default: first of `--metadata-cols`), `--target-cluster-size` (default 200, band [100, 300]), `--min-cluster-cells` (default 20), `--clustering-signal` (`wnn`/`rna`/`atac`; auto-selected from `sc_modality`), `--transferred-label-col` (default `predicted.id`), `--peak-caller` (`macs2`/`macs3`), `--rna-only`/`--atac-only`, `--dry-run`, `--yes`.

**Pipeline (orchestrated by `pseudobulk.py`):**
1. `detect-dataset-type` gate — refuses bulk; for sc-undetermined refuses with coembed-construct pointer; for separate-assay points at coembed-construct (run that first to produce a shared-space coembed.rds, then call this skill with `--use-existing-clusters`). The legacy ATAC-LSI clustering path with `--transferred-label-col` is still available for users who've integrated upstream themselves.
2. Dependency check (Rscript + MACS2 on PATH; degrades to warning under `--dry-run`).
3. `load_and_cluster.R` — loads object (extension dispatch: `.rds` direct, `.h5ad` via SeuratDisk *iff* the user installed it manually — no longer auto-installed, see SC notes in dependencies.yml — `.h5mu` via MuDataSeurat), runs Seurat/Signac WNN clustering (or RNA-only PCA / ATAC-only LSI), scale-aware resolution binary search (seed `r0 = 0.15·log2(N/1000)·sqrt(target_n_clusters/20)`, clamped to [0.05, 3.0], up to 8 iterations targeting mean cluster size in [100, 300]), drops <20-cell clusters and sub-groups, writes `clusters.csv`, `resolution_trace.json`, `groups_plan.json`, and per-group barcode files.
4. `aggregate_rna.R` — sums raw counts per `(cluster × metadata_col × value)` group, writes 2-column GeneQuant TSV (no header) per group.
5. Per-cluster `Signac::CallPeaks` (driven from R via `call_peaks.R`): subsets the object per group, restricts to ATAC-side cells (using a `meta.data$assay == 'ATAC'` flag if present, else `Cells(obj[[atac_assay]])`), reconciles barcodes against `fragments.tsv.gz` via suffix-stripping (≥95% overlap required, fail-loud otherwise), attaches a per-group `CreateFragmentObject`, and shells out to MACS2/MACS3 through Signac. Emits `<group>_peaks.rds` (Signac object) and `<group>.narrowPeak` (ENCODE 10-col).
6. Emits `manifest.tsv` with the exact columns `build-taiji-input --samples` expects: `name, cohort, group, rna_seq, atac_seq, genome`. Zero-glue handoff.

**Outputs:**
```
output_dir/
├── clusters.csv                    per-cell barcode + cluster + metadata
├── resolution_trace.json           audit log of the resolution binary search
├── groups_plan.json                per-(cluster × meta_col × value) plan
├── per_cluster_barcodes/<name>.txt one barcode per line (legacy from the pre-Signac MACS pipe; kept for downstream debugging / manual reruns)
├── rna/<name>.tsv                  GeneQuant TSV per group
├── atac/<name>.narrowPeak          per-cluster Signac::CallPeaks peaks (ENCODE 10-col)
├── atac/<name>_peaks.rds           per-cluster Signac peak GRanges (saveRDS)
├── qc/umap.png                     QC UMAP panels (clusters / assay / metadata cols); skipped via --no-plot
├── qc/umap_coords.csv              2D embedding + cluster labels for re-plotting outside R
└── manifest.tsv                    feeds build-taiji-input directly
```

**Status:** Python orchestrator + gates verified across detect-dataset-type's 11 fixtures (bulk redirect, sc-undetermined refusal, multiome auto-WNN, separate-assay auto-ATAC, dep-check error path, --dry-run graceful degradation). R scripts syntax-checked against Seurat v5 / Signac ≥1.12. **R is available** at `/usr/local/bin/Rscript` (4.2.0) in the sandbox, but the Seurat/Signac/EnsDb stack isn't pre-installed — install via `bin/install.sh --profile sc` on the machine you're actually running on. Real-data integration testing happens on the target machine (Mac or SLURM); the sandbox can probe `.rds` structure via base R `readRDS` to inspect class slots / metadata cols / barcode patterns without the full SC stack — useful for pre-flight diagnosis even when end-to-end runs aren't possible here.

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

**Purpose:** Idempotently stage the reference data Taiji needs into `<output_dir>/dependencies_data/`. FASTA + GTF are downloaded from GENCODE via the EMBL-EBI mirror; the MEME motif file is **CIS-BP** (vendored at `cisbp_database/`, no network — `cisbp_human_2.meme` for human, `cisBP_mouse.meme` for mouse). Reads a single `reference_manifest.yml` declaring URLs, target paths, gunzip flags, and approximate uncompressed sizes per genome. **HOCOMOCO is not supported** — locking the project to a single motif source keeps PageRank values comparable across runs.

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
├── GeneRanks.tsv               # headline output: per-TF PageRank score (~1100 rows for cisbp_human_2.meme on hg38)
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

**Executor auto-detection:** `--executor auto` (default) picks `slurm` when `sbatch` is on PATH, else `local`. SLURM mode materializes `taiji_array.sbatch` from `skills/taiji-runner/templates/taiji_array.sbatch.template`, submits with `sbatch -a 1-N%MAX` (default MAX=10 — matches the user's "10 max parallel" preference), polls `squeue` until the array completes, then validates outputs the same way as local mode (count `GeneRanks.tsv` per sample). Async mode via `--no-wait` returns the job ID immediately. Local mode stays sequential by default (`--parallel 1`) for memory safety; SLURM mode defaults to `--parallel 10`.

## Skill composition (end-to-end pipeline)

```
    fetch-references --genome hg38 --output dependencies_data/  (one-time per site)
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
            │    (.h5mu / cellranger-arc)   (one .rds RNA + one .rds ATAC)
            │           │                             │
            │           │                  coembed-construct ──┐
            │           │                  (Stuart/Signac      │
            │           │                   integrate_atac:    │
            │           │                   anchors → impute   │
            │           │                   RNA → merge → joint│
            │           │                   PCA/UMAP → cluster)│
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
                  (`taiji run --config <sample>_config.yml`)      │
                              │                                   │
                              ▼                                   │
              workflow-log.finalize(status) ──────────────────────┘
```

Each step gates on the previous one and produces inputs in the next step's expected format with **no glue code in between**. The `manifest.tsv` from `pseudobulk-construct` is a drop-in for `build-taiji-input --samples`.

## External dependencies (not vendored)

Trimmed to reflect that all inputs are pre-aligned/pre-quantified — no BWA, STAR, RSEM, samtools, Picard, or Singularity needed. Per-skill dep manifests live at `skills/*/dependencies.yml`; `bash bin/doctor.sh` verifies them.

**Always required:**
- Python ≥ 3.10
- Pinned Python deps in `environment.yml`: click, pydantic, pyyaml, jinja2, rich, pandas, openpyxl
- Rscript ≥ 4.2 (only for `pseudobulk-construct`)
- R packages (bioconda): Seurat ≥ 5.0, Signac ≥ 1.12, Matrix, GenomicRanges, dplyr, optparse, jsonlite
- R packages (GitHub-only, installed by `bin/postinstall.R`): MuDataSeurat (for `.h5mu`). SeuratDisk (for `.h5ad`) was removed from the auto-install — it's unmaintained upstream and brittle to build on cluster envs. `.h5ad` ingestion now requires either a manual `remotes::install_github("mojaveazure/seurat-disk")` or upstream conversion to `.rds`.
- MACS3 ≥ 3.0 (or MACS2 — `pseudobulk-construct` auto-detects whichever is on PATH)
- gunzip + awk (standard Unix)

**Taiji binary** — installed by `bin/install-taiji.sh --system <centos|ubuntu|macos>`:
- CentOS / RHEL: `taiji-CentOS-x86_64` (auto-fetched)
- Ubuntu / Debian: `taiji-Ubuntu-x86_64` (auto-fetched)
- macOS: `taiji-macOS-XX-XX` (suffix varies by macOS version; install-taiji.sh prints instructions)
- Singularity: explicitly out of scope (no official .sif published, can't install on this site's SLURM)

For Eunice's UCSD setup, the binary additionally lives at `/stg3/data1/eunice/Projects/nextflow/nf_taiji/binaries/` and may be symlinked from there if preferred.

## Environment management

Composable profiles let users install only what their dataset needs. Bulk-only users save ~3.5 GB and ~20 min by skipping the R + Seurat/Signac stack.

```bash
bash bin/install.sh --system <macos|centos|ubuntu>            # base only (default)
bash bin/install.sh --system macos --profile sc               # base + sc
bash bin/install.sh --system macos --profile full             # base + sc + dev

bash bin/auto-install.sh --data-dir <path> --system macos     # detect-driven
bash bin/auto-install.sh --data-dir <path> --system macos --dry-run  # plan only

micromamba activate taiji-agent
bash bin/doctor.sh --profile base                             # filter by profile
```

### Profiles

| Profile | What's in it | Enables | Disk | Time |
|---------|--------------|---------|------|------|
| `base` (default) | Python + click/pydantic/pyyaml + pandas + openpyxl + macs3 + local pkg | detect-dataset-type, build-taiji-input, fetch-references, taiji-runner, workflow-log | ~500 MB | ~5 min |
| `sc` (additive) | r-base + r-seurat + r-signac + Bioconductor (GenomicRanges + EnsDb) + r-remotes + (postinstall) MuDataSeurat | pseudobulk-construct, coembed-construct | +3-4 GB | +15-30 min |
| `dev` (orthogonal) | pytest + pytest-cov + ruff + mypy + ipython | author tooling | +500 MB | +3 min |
| `full` | base + sc + dev | everything | ~5 GB | ~25 min |

`sc` is **additive on top of base**, not standalone — you can layer profiles incrementally. A user who installed `base` and later needs SC can run `bash bin/install.sh --profile sc` and only the R packages get added; no full reinstall.

### Per-skill profile membership

Each skill declares its profile in `skills/<name>/dependencies.yml`:

```yaml
profile: base    # detect-dataset-type, build-taiji-input, fetch-references,
                 # taiji-runner, workflow-log
profile: sc      # pseudobulk-construct (the only SC-using skill)
```

`bin/doctor.sh --profile <name>` filters the verification table to skills that actually require the requested profile.

### Reproducibility tiers (full details in `docs/environment.md`)

1. **Tier 1 — micromamba + environment.<profile>.yml.** One command, ~5-25 min depending on profile.
2. **Tier 2 — conda-lock.** Generate `conda-lock.linux-64.yml` once, commit it, then `bash bin/install.sh --use-lockfile` produces byte-identical envs across machines.
3. **Tier 3 — conda-pack tarball.** After Tier 2 works on one machine, pack to a portable `.tar.gz`. Restore on any compatible Linux node in <2 min, no network. Closest equivalent to a Docker image without needing a container runtime.

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

## Single-cell preprocessing conventions (project-wide)

These apply to every `.rds` / `.h5ad` / `.h5mu` that flows into `coembed-construct` or `pseudobulk-construct`. The skills detect violations at startup and warn, but harmonizing upstream is strictly better than relying on the runtime warning — it eliminates the failure class entirely instead of just surfacing it.

- **Lowercase every metadata-value string before saving.** This is the single highest-value cleanup. The D7 dataset shipped with `tissue ∈ {spleen, siIEL}` on RNA and `tissue ∈ {Spleen, siIEL}` on ATAC; after `merge()` the column has THREE distinct values and a downstream `--metadata-cols tissue` cross-product creates separate groups for `spleen` vs `Spleen`, halving the spleen-cohort effective N with neither group flagged as wrong. The `coembed-construct` skill catches this with a `WARN: CASE-ONLY mismatch` message — but the user-side fix is still required to get a clean run. Idiomatic R, applied to BOTH per-modality objects before they enter the pipeline:
  ```r
  for (col in c("tissue", "genotype", "sample_id", "donor", "batch")) {
    if (col %in% colnames(obj@meta.data)) {
      obj[[col]][[1]] <- tolower(as.character(obj[[col]][[1]]))
    }
  }
  saveRDS(obj, "scRNA_clean.rds")
  ```
  This also collapses `WT`/`wt`, `Day7`/`day7`, `PBMC`/`pbmc`, etc. Free-text columns you don't pass through `--metadata-cols` are unaffected.
- **Use the same metadata column NAME for the same biological concept across modalities.** The skill validators only see what you ask them to compare. `sample_id` on RNA + `library_id` on ATAC for the same biological samples is invisible to `validate_metadata_consistency`. Pick a canonical name (`sample`, `donor`, `batch`) and use it on both.
- **Rebuild stale ChromatinAssay fragment paths before shipping `.rds` across machines.** `Fragments(atac)` embeds an absolute path at object-build time; on a destination machine the path resolves nowhere. The skill accepts `--fragments` as a runtime override, but the cleaner fix is to rebuild the Fragment handle on the destination machine before saving. See the `coembed-construct` SKILL.md "Recommended upstream preprocessing conventions" for the snippet.
- **Match the Seurat version of the input `.rds` to the env that consumes it.** Saving in v5 then loading under v4 silently produces an empty `Assays()` and crashes the pipeline ~10 min in. The skill refuses up-front under v4-with-v5-input, but the user-side fix (upgrade env or downgrade input via `as(obj[['RNA']], 'Assay')`) is preferable.

## Memory system

Persistent memory at `/sessions/quirky-festive-mayer/mnt/.auto-memory/`. Key files:

- `MEMORY.md` — index, loaded automatically.
- `user_role.md` — Eunice's role and preferences.
- `project_taiji_agent.md` — high-level project decisions.
- `project_detect_dataset_type.md` — skill-specific design notes.
- `project_pseudobulk_construct.md` — skill-specific design notes.

Read these before making non-trivial decisions about extending the project.

## Sandbox model (workspace-bound execution)

Every command the agent invokes runs through `bin/sandbox-run.sh`, a layered wrapper that ensures all writes stay inside the workspace folder. Three tiers, applied in order:

1. **Soft (always)**: cwd must be inside `REPO_ROOT`; `TMPDIR` is pinned to `<REPO_ROOT>/tmp/sandbox-<pid>`; child sees `$REPO_ROOT` for use in path scoping. Refuses to run if cwd is outside. Catches accidents (typos, agent mistakes) without OS support.
2. **OS-level (auto when available)**: `bwrap` on Linux or `sandbox-exec` on macOS. Filesystem writes are restricted to `REPO_ROOT` + workspace TMPDIR. Network is blocked unless `--allow-net` is passed. Engages automatically when the tool is present.
3. **Strict mode (opt-in)**: `--strict-sandbox` flag (or `TAIJI_STRICT_SANDBOX=1`) refuses to run unless an OS-level sandbox tool is present. Use on shared infrastructure where soft enforcement isn't enough.

Practical implications:

- `bin/run-taiji.sh` already wraps every Taiji invocation with `sandbox-run.sh`. Soft mode is the default; strict mode by env var or flag.
- `taiji_config.template.yml` uses `tmp_dir: ${REPO_ROOT}/tmp/<run_name>` instead of `/tmp/...`, so even Taiji's scratch I/O stays in the workspace and gets cleaned up alongside the run.
- The Taiji binary itself is x86_64 and the agent always shells out to it through the wrapper, so binary execution is workspace-bounded by inheritance.
- `--no-sandbox` exists for debugging; otherwise the wrapper is on by default.

The threat model this addresses is "agent makes a mistake or a buggy script writes outside the workspace," not "malicious code escapes the sandbox." For the latter, layer Docker/Apptainer on top — out of scope for this repo's design (no container runtime assumed on Mac/SLURM).

## Architectural constraints (load-bearing)

- **The Taiji binary cannot be executed from this sandbox.** Sandbox is aarch64; binaries are x86_64 ELF (Linux variants) or Mach-O (macOS). No qemu-user / box64 / fex-emu emulator is installed, and `sudo apt install` is blocked by `no_new_privileges`. The kernel rejects with `Exec format error`. Practical consequence: every "run Taiji" step is either Mac-side (Eunice's machine) or SLURM-side. The agent's role from the sandbox is strictly **prepare inputs + validate outputs**.
- **Workspace folder is shared between sandbox and Mac.** Files dropped into `<workspace>/data/`, `<workspace>/dependencies_data/`, `<workspace>/runs/<name>/Output/` from either side become immediately visible to the other. Use this for handoffs: Mac runs Taiji → outputs land in workspace → sandbox can validate / summarize / log without re-running.
- **R 4.2.0 IS in the sandbox** (`/usr/local/bin/Rscript`), but Seurat / Signac / EnsDb / GenomicRanges / optparse aren't pre-installed — Bioconductor compile-from-source takes ~30 min. MACS3 is also not in the sandbox by default. **What this means in practice:** you can use base R `readRDS` to inspect Seurat object slots / class / metadata columns / barcode formats WITHOUT installing the full SC stack — extremely useful for pre-flight diagnosis (catching issues like multi-sample barcode suffixes, pre-existing reductions, or unexpected assay layouts before submitting a long SLURM run). Full coembed/pseudobulk pipelines still belong on Mac/SLURM where the SC stack is installed.

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
