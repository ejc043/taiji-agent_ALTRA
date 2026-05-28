# pseudobulk-altra

Faithfully re-implements the pseudobulk construction block from `Taiji_ALTRA/scripts/prepare_taiji_input.r` (lines 183–356). Converts co-cluster assignments from a JointCCA object into per-cluster RNA GeneQuant TSVs and raw ATAC fragment BED.gz files — the exact format Taiji's `PairedEnd` ATAC pipeline expects.

Runs inside `singularity exec archr_1.0.3.sif` (for R) with data.table and Matrix.

## Inputs

| Argument | Type | Description |
|---|---|---|
| `--joint-cca` | path | `Save-Block1-JointCCA-cluster.rds` from archr-preprocess |
| `--labeled-rds` | path | Seurat RNA object with raw counts (labeled_v2.rds from archr-preprocess) |
| `--fragments` | path | Original fragments.tsv.gz |
| `--proj-name` | string | Sample prefix used in output filenames (e.g. `GSM8554115`) |
| `--output-dir` | path | Where pseudobulk files land |
| `--cluster-meta-dir` | path | Where cluster_identity.csv lands |
| `--min-cell` | int | Minimum cells per cluster for both RNA and ATAC (default: 20) |
| `--min-peak` | int | Minimum fragment rows per cluster for ATAC (default: 200) |

## Outputs

| File | Description |
|---|---|
| `<output-dir>/<proj>_C{i}_rna.tsv` | Two-column, no-header: `geneID<TAB>count` (raw counts summed across cluster) |
| `<output-dir>/<proj>_C{i}_atac.bed.gz` | Four-column BED: chr, start, end, barcode (raw fragments filtered to cluster cells, gzipped) |
| `<cluster-meta-dir>/<proj>_cluster_identity.csv` | Per-cluster: preClust (dominant cell type), purity, cells, RNA_cells, ATAC_cells |

Only clusters passing **both** `RNA_cells >= min-cell` AND `ATAC_cells >= min-cell` should be used as Taiji inputs. Clusters with only RNA TSV (no ATAC BED) or only ATAC BED (no RNA TSV) are excluded.

## Filter thresholds

| Filter | Default | Applied in |
|---|---|---|
| `min.cell = 20` | RNA cells ≥ 20 in cluster | `getCounts()` |
| `min.cell = 20` | ATAC cells ≥ 20 in cluster | `getPeaks()` |
| `min.peak = 200` | Fragment rows ≥ 200 | `getPeaks()` |

## Critical: cluster_meta alignment bug in original script

The original `prepare_taiji_input.r` computes `RNA_cells` and `ATAC_cells` via `lapply(unique(itg$Co_clusters), ...)` which returns values in the **order cells appear in `itg`**, then `cbind`s them against `preClust`/`cluster`/`purity`/`cells` which are in **sorted numeric cluster order** (from `table()`). This misaligns the RNA/ATAC cell counts with their cluster labels.

**Our fix**: iterate using `rownames(cM)` (the sorted cluster IDs) for both `RNA_cells` and `ATAC_cells`:

```r
RNA_cells  <- sapply(rownames(cM), function(x)
  sum(itg$Co_clusters == as.integer(x) & itg$Assay == "RNA"))
ATAC_cells <- sapply(rownames(cM), function(x)
  sum(itg$Co_clusters == as.integer(x) & itg$Assay == "ATAC"))
```

The cluster_meta CSV is for documentation only — the actual filtering happens inside `getCounts()` and `getPeaks()` which use the correct cell counts directly.

## `getPeaks` naming — raw fragments, not called peaks

Despite the name `getPeaks`, the function extracts **raw fragment rows** from the fragments.tsv.gz for each cluster's ATAC cells. It does **not** call MACS2/3. Taiji handles peak calling internally when it receives ATAC input with `tags: PairedEnd`.

## clusKNN ↔ FindClusters equivalence

`scDataviz::clusKNN()` is a one-line wrapper:

```r
clusKNN <- function(data, ...) {
  nn <- FindNeighbors(data, ...)          # k.param=20 (default)
  cl <- FindClusters(nn$snn, resolution=0.8, ...)
  as.integer(cl[,1])
}
```

Our implementation calls `FindNeighbors` and `FindClusters` identically with `set.seed(10)` immediately before, matching the original.

## Reference script

`data/ALTRA/pseudobulk_p1.R` — production script for Patient 1 (PB00072).
`data/ALTRA/cluster_meta/GSM8554115_cluster_identity.csv` — example output.
