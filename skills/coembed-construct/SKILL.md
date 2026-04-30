---
name: coembed-construct
description: Co-embed a separate-assay scRNA-seq + scATAC-seq pair into a shared latent space following the Stuart 2019 / Signac vignette workflow. Imputes RNA expression into ATAC cells via gene-activity-anchored TransferData, merges both objects, runs joint PCA + UMAP, and clusters on the shared space — producing a single Seurat .rds with `assay` (RNA/ATAC origin) and `seurat_clusters` (de novo, on the shared space) ready for pseudobulk-construct. Use this skill whenever the user has two separate Seurat objects (one RNA, one ATAC) for the same biological system and wants them in one coembed object — phrases like "co-embed RNA and ATAC", "integrate scRNA + scATAC", "WNN-style integration but for separate-assay", "Stuart vignette workflow", "joint embedding", or "shared UMAP across modalities".
---

# Co-embed scRNA-seq + scATAC-seq

## What this skill does

Implements the [Stuart et al. 2019 / Signac integrate_atac vignette](https://stuartlab.org/signac/articles/integrate_atac) workflow:

1. **Modality preprocessing** — RNA: NormalizeData → FindVariableFeatures → ScaleData → RunPCA → RunUMAP. ATAC: RunTFIDF → FindTopFeatures → RunSVD → RunUMAP.
2. **GeneActivity** — Signac computes per-gene ATAC accessibility (2 kb upstream + gene body) → adds an `ACTIVITY` assay on the ATAC object.
3. **FindTransferAnchors** — RNA reference, ATAC query (using ACTIVITY as query.assay), CCA reduction.
4. **TransferData on RNA expression** — imputes a full RNA expression matrix for every ATAC cell using the gene-activity anchors. Adds the imputed matrix as the `RNA` assay on the ATAC object.
5. **Merge** — `merge(rna, atac)` produces a single Seurat object with cells from both modalities.
6. **Joint PCA + UMAP** — ScaleData / RunPCA / RunUMAP on the merged object using only RNA features (real for RNA-origin cells, imputed for ATAC-origin cells).
7. **Cluster on the shared space** — FindNeighbors + FindClusters with the same scale-aware resolution binary search the pseudobulk-construct skill uses (target ~200 cells/cluster, range [100, 300], drops <20-cell clusters).
8. **Tag origin** — adds `meta.data$assay` ∈ {"RNA", "ATAC"} so downstream skills can filter to RNA-origin cells for counts aggregation and ATAC-origin cells for peak calling.

## What this skill does NOT do

- **Cell-type label transfer** (`TransferData(refdata = seurat_annotations)` from the vignette). Cluster IDs come from de novo clustering on the shared space; we don't project cell-type labels from a reference. If you want named cell types, do label transfer yourself upstream and they'll be carried through as a metadata column.
- **Sample integration / batch correction**. If you have multiple samples per modality, integrate them into one RNA object and one ATAC object first (Harmony / RPCA / scVI / etc.); coembed-construct then combines those two integrated objects.
- **Peak calling, RNA aggregation, manifest generation**. That's `pseudobulk-construct`'s job — coembed-construct's output feeds it directly.

## Inputs

| Flag | Purpose |
|------|---------|
| `--rna` | Path to scRNA-seq Seurat `.rds`. Must have an `RNA` assay with raw counts. Required. |
| `--atac` | Path to scATAC-seq Seurat `.rds`. Must have a `ChromatinAssay` (named `ATAC` or `peaks`) with a fragments handle. Required. |
| `--output` | Where to write the coembed `.rds`. Required. |
| `--genome` | hg38 / hg19 / mm10. Drives the EnsDb annotation pulled in for `GeneActivity`. Required. |
| `--target-cluster-size` | Target mean cluster size for the resolution binary search. Default: `200` (range `[100, 300]`). |
| `--min-cluster-cells` | Drop clusters with fewer than this many cells. Default: `20`. |
| `--metadata-cols` | Comma-separated metadata cols to preserve in the output. Skill checks they exist on **both** input objects (or warns if not) — values from the original objects are merged through verbatim. |
| `--n-pcs` | Number of PCs to use for joint UMAP / FindNeighbors. Default: `30`. |
| `--lsi-skip-first` | Skip the first N LSI components on the ATAC side (typically depth-correlated). Default: `1`. |
| `--no-cluster` | Stop after joint UMAP; don't cluster. Useful if the user wants to cluster manually. |
| `--no-plot` | Skip the QC UMAP rendering. |
| `--resolution` | Force a single resolution instead of binary search. Default: unset (binary search runs). |
| `--reuse-rna-reductions` | If the RNA object already has `pca` + `umap` reductions plus VFs and a non-empty `data` slot, skip the standard preprocessing pipeline (Norm/VF/Scale/PCA/UMAP). Saves 5–10 min on 60k+ cells. Falls back to full preprocessing if any prerequisite is missing. |
| `--reuse-atac-reductions` | If the ATAC object already has `lsi` + `umap.atac` (or `umap`) reductions, skip TF-IDF/SVD/UMAP. Saves 3–5 min on 20k+ cells. |
| `--strict-metadata` | Fail-loud if any `--metadata-cols` values differ between RNA and ATAC inputs (e.g. tissue=`spleen` on RNA, tissue=`Spleen` on ATAC). Default: warn but continue. Use this for production runs where silent group-splitting downstream would corrupt your stratification. |
| `--abort-on-memory-risk` | Refuse to run when the macOS pre-flight estimates peak RAM > available RAM + 0.5×free_swap (the regime where jetsam-OOM during CCA is likely). Default: warn but continue. Use on shared Macs where a 20-min thrash followed by a silent kill is worse than a fail-fast at startup. |
| `--fragments` | Path to `fragments.tsv.gz`. Required when the path embedded in the input ATAC `.rds` is stale — common when the object was built on HPC and shipped to a laptop, since `Fragments(atac)` carries an absolute path that won't resolve on the destination. Skill rebuilds the Fragment handle in place. Unused if the embedded path resolves. |

## Outputs

```
output_dir/
├── coembed.rds                # Seurat object: cells from both RNA + ATAC, shared PCA/UMAP, seurat_clusters from the shared space, meta.data$assay tags origin
├── coembed_summary.json       # chosen resolution, n_clusters, n_cells per (cluster × assay), genes_used count, anchor count
└── qc/
    ├── umap.png               # joint UMAP, ALWAYS rendered (unless --no-plot). See panel contract below.
    ├── umap_pre_merge.png     # Stuart-vignette-style side-by-side: RNA's own UMAP next to ATAC's own UMAP, BEFORE integration. The "did integration converge well?" reference plot.
    └── umap_coords.csv        # 2D embedding coords + cluster + assay + every --metadata-cols value (for re-plotting in matplotlib/plotly outside R)
```

### `qc/umap.png` panel contract

Generated by default (skip with `--no-plot`). Panels rendered, in order:

1. **`seurat_clusters`** — de novo clusters on the shared joint PCA. Title shows `(res=<chosen>, n=<n_clusters>)`. Skipped if `--no-cluster`.
2. **`assay`** — RNA-origin vs ATAC-origin cells, colored with explicit qualitative colors (RNA = `#1f77b4` blue, ATAC = `#ff7f0e` orange) so the distinction is unmistakable on grayscale prints or for colorblind viewers. Title shows cell counts: `assay (RNA: <n>  /  ATAC: <n>)`. **This is the integration-QC panel — well-integrated data shows both colors mixed across clusters; poorly-integrated data shows clusters segregating by assay.**
3. **One panel per `--metadata-cols` entry**, in the order the user passed them. The skill validates every requested column against the merged object's `meta.data` at plot time and emits a `WARN: --metadata-cols not present on merged object (skipped panels): <list>` if any are missing.

Layout: 2-column grid via `patchwork::wrap_plots(ncol = 2)`. So with `--metadata-cols genotype,tissue` you get a 2×2 grid: clusters / assay / genotype / tissue.

### `qc/umap_pre_merge.png` panel contract

Two panels, side-by-side:
- **RNA only** — the RNA object's standalone UMAP (computed at stage 1 of the pipeline before any integration), colored solid blue (`#1f77b4`), titled `RNA only (n=<n_cells>)`.
- **ATAC only** — the ATAC object's standalone UMAP (computed at stage 2 from LSI dims 2:30), colored solid orange (`#ff7f0e`), titled `ATAC only (n=<n_cells>)`.

This is the Stuart vignette's [first figure](https://stuartlab.org/signac/articles/integrate_atac) — you want it as a sanity check before trusting the joint UMAP. If either modality looks structureless on its own, the joint embedding inherits that problem.

### Where to put `--output` (convention)

**Always under `runs/<run_name>/coembed/coembed.rds`, not inside the input data directory.** The skill emits a `WARN` on stderr if it detects the output going under the RNA/ATAC parent dir. Reasons:

- Cleanup is `rm -rf runs/<run_name>/`. With outputs in `data/<dataset>/coembed/`, you'd have to manually scrub.
- Re-runs with different parameters (resolution, metadata cols) need separate output dirs. Pinning to a run dir gives you that automatically.
- A coembed.rds dropped into the input data dir confuses the next `detect-dataset-type` scan — the gate sees the original two `.rds` plus the coembed and may produce ambiguous classifications. Keeping the input dir immutable avoids this entirely.

```bash
# Right:
python skills/coembed-construct/scripts/coembed.py \
  --rna  data/<dataset>/rna.rds \
  --atac data/<dataset>/atac.rds \
  --output runs/<run_name>/coembed/coembed.rds \
  --genome mm10

# Wrong (skill warns):
  --output data/<dataset>/coembed/coembed.rds   # <- inside input dir
```

The output `coembed.rds` is the input to `pseudobulk-construct --input <path>/coembed.rds --use-existing-clusters` (since clustering already happened here).

## Real-data gotchas and how the skill handles them

These come up in practice on multi-sample 10x datasets. Each is now logged at startup so the user sees the situation before the long-running steps fire.

| Gotcha | What it looks like | Skill behavior |
|--------|-------------------|----------------|
| **Multi-sample-aggregated barcodes** with `-1`/`-2`/.../`-N` suffixes (cellranger-aggr / 10x Multi output) | `AAACCAAAGAACGGCA-1`, `AAACGAAAGAACGACC-2` etc. The fragments file carries all sample suffixes; each Seurat object holds the subset that passed QC. | Identity-match works as-is; my barcode reconciliation tries identity FIRST and only suffix-strips if that fails. No action needed — just don't be surprised when `gunzip -c fragments.tsv.gz | head` shows multiple suffixes. |
| **Pre-existing reductions** on inputs (`pca, umap` on RNA / `lsi, umap` on ATAC) | The user has already run their own per-modality clustering and embedding before handing objects off. | Skill RE-COMPUTES from raw counts by default (reproducibility). Pass `--reuse-rna-reductions` and/or `--reuse-atac-reductions` to skip — saves ~10 min on a 60k+21k pair. Logged at startup: `RNA input: 60000 cells, reductions={pca,umap}`. |
| **Pre-existing `seurat_clusters`** column on either input | The user already clustered each modality. After joint clustering, the original labels would be silently overwritten. | Skill ALWAYS preserves them under `rna_input_clusters` and `atac_input_clusters` BEFORE merging, so they survive into the coembed object's `meta.data`. You can compare original-per-modality vs. de-novo-joint clusters downstream. |
| **Pre-existing `RNA` assay on the ATAC object** | E.g. from a prior `integrate_atac` run, or because the ATAC object was sliced from a multiome. | Skill OVERWRITES it with freshly imputed values from `TransferData`. Logged loudly: `NOTE: ATAC object already carries an 'RNA' assay. It will be REPLACED later`. If you want to keep the existing values, save them under a different assay name (e.g. `atac[["RNA_old"]] <- atac[["RNA"]]`) before invoking this skill. |
| **Large cell counts** (>50k merged) | Joint PCA/UMAP on 60k+ cells with ~2k variable features needs 10–20 GB peak RAM. | No code change — just an SBATCH `--mem` reality check. Recommend `--mem=128G` for >80k merged cells. |
| **Metadata cols missing on one input** | e.g. you pass `--metadata-cols genotype,batch_id` but `batch_id` only exists on the ATAC side. | After merge, the column is NA for cells from the side that lacked it. Plotting still works (NA shows as gray). The QC plot logs a `WARN` listing skipped cols. |
| **Metadata VALUE / CASING mismatch across inputs** | RNA has `tissue ∈ {spleen, siIEL}`, ATAC has `tissue ∈ {Spleen, siIEL}`. After merge, `tissue` has THREE distinct values; downstream `pseudobulk-construct --metadata-cols tissue` would create separate groups for `spleen` vs `Spleen`. | Skill validates `--metadata-cols` value sets across both inputs at startup. Logs a `WARN: CASE-ONLY mismatch in 'tissue'` (or `WARN: VALUE mismatch`) listing the offending values, with a copy-paste fix hint (e.g. `atac$tissue <- tolower(atac$tissue)`). Pass `--strict-metadata` to fail-loud instead of warning. |
| **Seurat v5 (`Assay5`) inputs in a Seurat v4 env** ⚠️ **HARDEST silent failure** | `assays$RNA` is class `Assay5` (Seurat v5) but the running Seurat is v4.x. `readRDS` succeeds, `class(obj)` is "Seurat", but `Assays(obj)` returns empty and `GetAssayData(obj, slot="counts")` errors with `'RNA' is not an assay` ~10 min into the pipeline. | Skill checks `packageVersion("Seurat")` at startup and the per-assay class on each input. Refuses with: `RNA assay 'RNA' is class Assay5 but Seurat <ver> is loaded. ... Either upgrade the env to Seurat >= 5, or downgrade the input via obj[['RNA']] <- as(obj[['RNA']], 'Assay'); saveRDS(...)`. Verified on Seurat 4.3.0 reading the D7 dataset's v5 RNA object. |
| **Seurat v5 (`Assay5`) inputs in a Seurat v5 env** | Counts live in `attributes(.)$layers$counts`, not the v3 slot. Raw `attributes()` probes look empty but `GetAssayData()` dispatches correctly. | No code change needed; `GetAssayData()` / `NormalizeData()` / `RunPCA()` all dispatch across v3/v5 in Seurat >= 5. The `--reuse-rna-reductions` check uses `GetAssayData(slot="data")` which works on v5 layers. |
| **Mixed-version RNA/ATAC** | RNA is Seurat v5 (Assay5), ATAC is Seurat v3 (Assay / ChromatinAssay). | `merge()` handles this in Seurat ≥ 5.0; the result follows the active class on the LHS object. The skill loads RNA first so a v5 RNA produces a v5 coembed object. If `merge` errors on your specific input, downgrade RNA via `JoinLayers(rna)` upstream. |
| **Stale fragments path inside ChromatinAssay** | The `.rds` was built on HPC and shipped to a laptop. `Fragments(atac)` carries the absolute path from the build machine (e.g. `/home/whw006/work/.../fragments.tsv.gz`); on the destination machine `file.exists()` returns FALSE. `readRDS` succeeds; `GeneActivity()` then errors ~10 min later with an opaque tabix/seek error. | Skill checks each Fragment object's `path` against the local filesystem at startup and logs `paths exist=TRUE,FALSE,...`. If any are missing AND `--fragments` is provided, the skill rebuilds a single Fragment object in place (preserving the existing barcode-cell mapping). If any are missing AND `--fragments` is NOT provided, the skill refuses loudly with a copy-paste fix. |
| **`GetAssayData(slot=...)` defunct in SeuratObject ≥ 5.0** | SeuratObject 5.0 deprecated `slot=` and 5.0+ made it a hard error. Code that worked under Seurat v4 (`GetAssayData(rna, slot = "counts")`) errors immediately under v5 with a lifecycle stop. | Skill uses `layer = ...` everywhere. Stays compatible with the v3/v5 dispatch. |
| **macOS R `R_MAX_VSIZE` cap (~16-18 GB) below what FindTransferAnchors needs** | On Apple Silicon Macs, R imposes a default vector-memory cap independent of physical RAM. CCA between a 60k RNA reference and a 21k ATAC query × 2k features needs ~10-20 GB working memory; R errors with `vector memory limit of 18.0 Gb reached, see mem.maxVSize()` mid-anchoring. Linux R has no equivalent cap. | Skill detects `Sys.info()$sysname == "Darwin"` at startup and bumps the cap to 128 GB (override via `TAIJI_R_MAX_VSIZE_GB` env var). On Macs with less than ~24 GB physical RAM, the kernel will swap the spike to NVMe SSD — slower but completes. On 16-18 GB Macs the run still finishes (~30 min) thanks to swap. SLURM nodes need no change. |
| **macOS jetsam OOM-kill mid-CCA (the *swap-budget-exhausted* failure)** | Bumping VSIZE prevents R from self-throttling, but doesn't stop the macOS kernel from killing R when physical RAM AND swap budget are both saturated. On the D7 dataset (60k RNA + 21k ATAC × 1639 features) running on an 18 GB Mac, FindTransferAnchors peaked at ~10 GB RSS, swap filled to 23.7/24.6 GB (97%), the process went into uninterruptible-sleep (state `U`) thrashing at ~3% CPU, and progress stalled. Whether the run *eventually* finishes or gets jetsam-killed depends on the moment-to-moment swap pressure — the kernel kills R via SIGKILL (or sometimes SIGTERM) with **NO R error message**; the user just sees `coembed.R exited -15` after 5+ wasted minutes. | Skill computes a pessimistic peak-RSS estimate `8 * n_features * (n_rna + n_atac) * 6 bytes` after both objects load (before GeneActivity), looks up `hw.memsize` + `vm.swapusage` via sysctl, and emits `INFO` / `WARN` / `HARD-WARN` based on `peak_estimate vs RAM + 0.5*free_swap`. Pass `--abort-on-memory-risk` to refuse instead of warning. The Python wrapper also recognizes signal-driven exits (-9 / -15 / 137 / 143) at the Rscript layer and prints a multi-line OOM diagnostic with concrete mitigations: run on SLURM, `subset(rna, downsample = 20000)`, lower `--n-pcs`, free swap. The estimate is captured in `coembed_summary.json` as `peak_ram_estimate_gb` for post-run audit. |
| **High `FindTransferAnchors` feature drop rate (silent quality issue)** | `FindTransferAnchors(features = VariableFeatures(rna), ...)` silently drops every feature absent from EITHER assay (the RNA reference OR the ACTIVITY assay). Signac's `GeneActivity` only computes coverage for genes with EnsDb annotations, so a baseline ~10-20% drop is normal. But >50% drop usually signals a real upstream problem: wrong `--genome` for the input data, EnsDb species mismatch (e.g. mm10 EnsDb on human RNA), or gene-name convention mismatch (Ensembl IDs vs symbols). Currently FindTransferAnchors emits a one-line stderr warning that's easy to miss inside Seurat's normal chatter. | Skill computes `intersect(VariableFeatures(rna), rownames(atac[["ACTIVITY"]]))` BEFORE FindTransferAnchors and logs `features for anchoring: N requested -> M shared (X% kept, Y dropped)`. Drop rate >50% triggers an explicit `WARN` with the most likely root causes. Counts are written to `coembed_summary.json` as `n_features_requested` / `n_features_used` / `n_features_dropped` for auditing across runs. On the D7 mm10 dataset, 1639/2000 (82%) features were retained — a healthy rate confirming the genome / annotation / gene-name pipeline is consistent. |
| **Louvain modularity saturation in resolution search** | On 60k+ merged cells with `--target-cluster-size 200`, Louvain plateaus around 75-100 clusters (~1000 cells/cluster mean) and refuses to split further regardless of resolution. The binary search bumps the resolution to the upper bracket (`hi=5.0`), gets the same partition, and if naively run, burns 8 iters of `FindClusters` (~60s each on 81k cells) producing identical work. | Search detects when `n_clusters` is unchanged across two consecutive iters at the same `hi` (or `lo`) bracket bound and breaks out early with `[coembed] resolution search hit Louvain modularity ceiling at r=...`. Saves ~5-7 min per run on the saturation case. The chosen resolution is still the saturating one; downstream stratification works fine on the resulting clusters. If you genuinely want finer granularity, sub-cluster the desired parent cluster with `FindSubCluster()` post hoc. |

## Pre-flight checks the skill enforces

- **RNA object** must have `RNA` assay with raw counts (`GetAssayData(rna, layer = "counts")` non-empty).
- **ATAC object** must have a ChromatinAssay (named `ATAC`, `peaks`, or any class-`ChromatinAssay` assay) with fragments.
- **Fragments-path resolves on the local machine** — every Fragment object's embedded `path` slot is `file.exists()`-checked. If any are missing the skill demands a `--fragments` override (or refuses).
- **Genome support** — fails loud on unsupported `--genome` values; the EnsDb dispatch table maps human/mouse only.
- **Metadata-cols sanity** — if a column is in `--metadata-cols` but missing from one input, prints a warning and continues (the merged object will have NA for the missing-side cells in that column).

## Recommended upstream preprocessing conventions

The skill validates a lot of failure modes at startup, but a few conventions on the user-side `.rds` build prevent the worst silent bugs from ever materializing in the first place. Apply these when you generate the per-modality objects (Seurat / Signac scripts that produce `*_scRNA.rds` / `*_scATAC.rds`):

- **Lowercase every metadata-value string before saving the `.rds`.** This is the single biggest unforced-error fix. The D7 dataset shipped with `tissue ∈ {spleen, siIEL}` on RNA but `{Spleen, siIEL}` on ATAC; after `merge()` the column has THREE distinct values (`spleen`, `Spleen`, `siIEL`), and a downstream `pseudobulk-construct --metadata-cols tissue` then creates separate groups for `spleen` and `Spleen` — silently halving the effective sample size of one tissue. The skill catches this at startup with a `WARN: CASE-ONLY mismatch` message, but harmonizing upstream is strictly better than relying on the warning. Idiomatic R:
  ```r
  for (col in c("tissue", "genotype", "sample_id", "donor", "batch")) {
    if (col %in% colnames(obj@meta.data)) {
      obj[[col]][[1]] <- tolower(as.character(obj[[col]][[1]]))
    }
  }
  saveRDS(obj, "scRNA_clean.rds")
  ```
  Apply the same loop to BOTH the RNA and ATAC objects before they ever enter `coembed-construct`. Free-text values that happen to match across modalities (e.g. `"WT"` vs `"wt"`, `"day7"` vs `"Day7"`, `"PBMC"` vs `"pbmc"`) all collapse cleanly under `tolower()` and won't split groups.
- **Use the SAME column NAME for the same biological concept across RNA + ATAC.** Real example from D7: RNA had `sample_id` ∈ {`siIEL_KO_D7`, …}; ATAC had `library_id` ∈ {`siIEL_KO_D7_CRISPR_18`, …}. Same concept, different column names AND different value formats — neither the metadata-value validator nor `--metadata-cols` can rescue this; they only see what you ask them to compare. Pick one canonical column name (e.g. `sample`) on the upstream pipeline and rename both modalities to match.
- **Never carry stale path strings across machines.** `Fragments(atac)` embeds the absolute fragments path at object-build time. When a `.rds` ships from HPC to a laptop, the path resolves nowhere. The skill detects this and accepts `--fragments` to rebuild in place — but if you'd rather avoid the runtime kludge, rebuild the Fragment object on the destination machine before saving:
  ```r
  fr <- CreateFragmentObject(path = "/local/path/fragments.tsv.gz",
                             cells = colnames(atac))
  Fragments(atac[["peaks"]]) <- NULL
  Fragments(atac[["peaks"]]) <- fr
  saveRDS(atac, "scATAC_clean.rds")
  ```
- **Don't ship a Seurat v5 `.rds` to a v4-pinned env.** The skill detects the `Assay5`-under-Seurat-4 trap at startup and refuses, but you save 30 min by either upgrading the destination env to Seurat ≥ 5.0 OR downgrading the input via `obj[['RNA']] <- as(obj[['RNA']], 'Assay'); saveRDS(...)` BEFORE shipping.

## When to run this skill vs. pseudobulk-construct directly

| Input | Skill |
|-------|-------|
| Two separate `.rds` (one RNA, one ATAC), same biological system | **coembed-construct** first → its output → pseudobulk-construct |
| One pre-coembedded `.rds` (you ran `merge` + clustering yourself, has `assay` + `seurat_clusters` columns) | pseudobulk-construct directly with `--use-existing-clusters` |
| One `.h5mu` (multiome — same cells in both modalities) | pseudobulk-construct directly with `--clustering-signal wnn` (no coembed needed; cells are already paired) |
| One `.rds` (RNA-only or ATAC-only) | pseudobulk-construct directly with `--clustering-signal rna` or `atac` |

## Why we don't auto-do label transfer

Cluster IDs from de novo clustering on the shared space are sufficient for the cross-product stratification (cluster × genotype × tissue) the downstream `pseudobulk-construct` does. Reference-derived cell-type labels (`TransferData(refdata = seurat_annotations)` from the vignette) are a separate decision — which reference atlas to anchor against affects which labels land on which cells, and bad anchors silently produce confidently-wrong labels (worst failure mode for SC analysis). If you want named labels, run `TransferData` yourself upstream and pass the resulting column through `--metadata-cols`.

## Composition with the rest of the agent

```
detect-dataset-type → separate-assay
                          │
                          ▼
   coembed-construct (--rna, --atac, --output)
                          │
                          ▼
              <output>/coembed.rds
                          │
                          ▼
   pseudobulk-construct (--input coembed.rds --use-existing-clusters --fragments ...)
                          │
                          ▼
                    manifest.tsv
                          │
                          ▼
   build-taiji-input → taiji-runner → GeneRanks.tsv per pseudobulk
```

The skill auto-attaches to `workflow-log` if a run is active; logs anchor count, gene-list size, chosen resolution, and per-cluster RNA/ATAC cell counts.

## Where to find more detail

- `references/coembed_strategy.md` — design notes, mapping to the Stuart vignette, why we use CCA / which dims, how the resolution search adapts to coembedded cell counts.
