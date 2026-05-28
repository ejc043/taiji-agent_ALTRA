# archr-preprocess

Faithfully re-implements the ArchR preprocessing block from `Taiji_ALTRA/scripts/prepare_taiji_input.r` for a single scATAC-seq sample. Produces an ArchR project with IterativeLSI, UMAP, predicted cell-type labels from Seurat reference transfer, and a JointCCA co-embedding ready for pseudobulk construction.

Runs inside `singularity exec archr_1.0.3.sif`.

## Inputs

| Argument | Type | Description |
|---|---|---|
| `--fragments` | path | fragments.tsv.gz (bgzipped, tabix-indexed) |
| `--proj-name` | string | Sample identifier — used as Arrow prefix and ArchR project name |
| `--output-dir` | path | Where to write the ArchR project directory |
| `--rna-h5` | path | HISE-format HDF5 RNA file for the matched RNA sample |
| `--pbmc-ref-rds` | path | pbmc_reference.rds (pre-converted from pbmc_multimodal.h5seurat, Platelet-filtered) |
| `--atac-barcodes` | path | TSV/CSV with a `barcodes` column; only these ATAC barcodes are kept (PassQC=TRUE, singlet=TRUE from backup_meta or equivalent) |
| `--r-libs` | path | Extra R library path (default: `/stg3/data1/eunice/R_libs/archr_extra`) |
| `--threads` | int | ArchR + Seurat threads (default: 5) |

## Outputs

| File | Description |
|---|---|
| `<output-dir>/` | ArchR project (Arrow, cellColData with predictedGroup_Un, Co_clusters) |
| `<output-dir>/RNAIntegration/GeneIntegrationMatrix/Save-Block1-JointCCA.rds` | Raw CCA embedding (all cells × CC dims + metadata) |
| `<output-dir>/RNAIntegration/GeneIntegrationMatrix/Save-Block1-JointCCA-cluster.rds` | CCA + Co_clusters column |
| `<output-dir>/archr_umap_done.flag` | Checkpoint: IterativeLSI→UMAP done |
| `<output-dir>/archr_intg_done.flag` | Checkpoint: gene integration + clusKNN done |
| `<proj-name>_labeled_v2.rds` | Seurat RNA object with `ref.spca` reduction and `predicted.celltype.l2` |

## Key parameters (verbatim from prepare_taiji_input.r)

| Parameter | Value | Stage |
|---|---|---|
| `minTSS` | 4 | createArrowFiles |
| `minFrags` | 1000 | createArrowFiles |
| `iterations` | 4 | addIterativeLSI |
| `varFeatures` | 25000 | addIterativeLSI |
| `sampleCellsPre` | 10000 | addIterativeLSI |
| `projectCellsPre` | FALSE | addIterativeLSI |
| `sampleCellsFinal` | 100000 | addIterativeLSI |
| `resolution` | 3 | addClusters |
| `normalization.method` | "SCT" | FindTransferAnchors |
| `dims` | 1:50 | FindTransferAnchors |
| `reference.reduction` | "spca" | FindTransferAnchors |
| `recompute.residuals` | FALSE | FindTransferAnchors |
| `new.reduction.name` | "ref.spca" | IntegrateEmbeddings |
| `sampleCellsATAC` | 50000 | addGeneIntegrationMatrix |
| `sampleCellsRNA` | 50000 | addGeneIntegrationMatrix |
| `groupRNA` | "predicted.celltype.l2" | addGeneIntegrationMatrix |
| `resolution` | 0.8 | FindClusters (clusKNN equivalent) |
| `set.seed` | 10 (before clusKNN) | clusKNN |

## Container requirements

```
singularity exec /stg3/data1/eunice/bin/containers/archr_1.0.3.sif Rscript run_archr_preprocess.R ...
```

Extra R libs needed (not in container): `TxDb.Hsapiens.UCSC.hg38.refGene`, `org.Hs.eg.db` — install to `--r-libs` path via `BiocManager::install()`.

## Deviations from original (unavoidable)

| Original | Our implementation | Reason |
|---|---|---|
| Pre-built HISE Arrow files | `createArrowFiles()` from fragments.tsv.gz | HISE cloud not available |
| `filterDoublets(filterRatio=0.5)` | Barcode whitelist from backup_meta (PassQC=TRUE & singlet=TRUE) | User-specified; equivalent quality gate |
| `LoadH5Seurat("pbmc_multimodal.h5seurat")` | `readRDS(pbmc_reference.rds)` | SeuratDisk not in container |
| `H5weaver::read_h5_seurat()` | Custom `read_h5_seurat_rhdf5()` | H5weaver not in container |
| `scDataviz::clusKNN()` | `FindNeighbors() + FindClusters(resolution=0.8)` | scDataviz not in container; defaults are identical |
| `addArchRGenome('hg38')` | `data(genomeAnnoHg38); GENOME_ANNO$genome <- "nullGenome"` | BSgenome.Hsapiens.UCSC.hg38 not in container |
| No seed before IterativeLSI | `set.seed(1)` added | Reproducibility |

## Gotchas

- **Seurat v5 Assay5 bug**: `options(Seurat.object.assay.version = "v3")` before `CreateSeuratObject` — the container's SeuratObject version has `Assay5@cells` as a character vector, causing `LayerData<-.Assay5` to fail on >1 cell.
- **Duplicate gene symbols**: HISE H5 files have non-unique gene names; apply `make.unique()` after reading.
- **`seqlevelsStyle` network call**: ArchR's gene annotation step tries to fetch chromInfo from UCSC. HPC nodes without internet need the `tryCatch` fallback in the script.
- **PROJ_NAME# prefix**: ArchR cellNames include the project prefix (`GSM8554115#barcode`). Strip with `sub(paste0("^", PROJ_NAME, "#"), "", cellNames)` before matching to external barcode lists.

## Reference script

`data/ALTRA/generate_p1_archr.R` — the production script used for Patient 1 (PB00072).
