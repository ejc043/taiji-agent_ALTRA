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

## Outputs

```
output_dir/
├── coembed.rds            # Seurat object: cells from both RNA + ATAC, shared PCA/UMAP, seurat_clusters from the shared space, meta.data$assay tags origin
├── coembed_summary.json   # chosen resolution, n_clusters, n_cells per (cluster × assay), genes_used count, anchor count
└── qc/
    ├── umap.png           # combined UMAP panels: clusters, assay (RNA vs ATAC), each --metadata-cols value
    └── umap_coords.csv    # 2D embedding coords + cluster + assay (re-plotting outside R)
```

The output `coembed.rds` is the input to `pseudobulk-construct --input coembed.rds --use-existing-clusters` (since clustering already happened here).

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
