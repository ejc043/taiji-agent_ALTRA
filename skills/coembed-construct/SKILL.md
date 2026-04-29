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
| `--genome` | hg38 / hg19 / mm10 / mm39. Drives the EnsDb annotation pulled in for `GeneActivity`. Required. |
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
| **Seurat v5 (`Assay5`) inputs with layer-based storage** | `assays$RNA` is class `Assay5` and counts/data live in `attributes(.)$layers$counts` rather than the v3 slot. `attributes(obj)$assays$RNA` looks empty in raw probes but `GetAssayData()` returns the data correctly. | The skill uses `GetAssayData()` / `NormalizeData()` / `RunPCA()` consistently — these all dispatch correctly across v3 / v5. No code change needed; just don't be alarmed when `str(obj)` shows empty slot/layer in old debug habits. The `--reuse-rna-reductions` check uses `GetAssayData(slot="data")` which works on v5. |
| **Mixed-version RNA/ATAC** | RNA is Seurat v5 (Assay5), ATAC is Seurat v3 (Assay / ChromatinAssay). | `merge()` handles this in Seurat ≥ 5.0; the result follows the active class on the LHS object. The skill loads RNA first so a v5 RNA produces a v5 coembed object. If `merge` errors on your specific input, downgrade RNA via `JoinLayers(rna)` upstream. |

## Pre-flight checks the skill enforces

- **RNA object** must have `RNA` assay with raw counts (`GetAssayData(rna, slot = "counts")` non-empty).
- **ATAC object** must have a ChromatinAssay (named `ATAC`, `peaks`, or any class-`ChromatinAssay` assay) with fragments.
- **Genome support** — fails loud on unsupported `--genome` values; the EnsDb dispatch table maps human/mouse only.
- **Metadata-cols sanity** — if a column is in `--metadata-cols` but missing from one input, prints a warning and continues (the merged object will have NA for the missing-side cells in that column).

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
