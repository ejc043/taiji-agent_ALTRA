#!/usr/bin/env Rscript
# aggregate_rna.R
#
# Sum raw RNA counts per (cluster x metadata_col x metadata_value) triple and
# write each as a GeneQuant TSV (gene_id<TAB>count) in the layout
# build-taiji-input expects.
#
# Why raw counts: bulk RNA-seq tools normalize internally and assume an
# integer count matrix. Summing normalized or log-transformed values gives
# meaningless "bulk" numbers.
#
# Reads:
#   --input       single-cell object (same one used by load_and_cluster.R)
#   --clusters    clusters.csv produced by load_and_cluster.R
#   --groups      groups_plan.json produced by load_and_cluster.R
#   --output-dir  where to write <group_name>.tsv files
#
# Writes:
#   <output-dir>/<group_name>.tsv   two columns, no header (GeneQuant format)

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(jsonlite)
})

option_list <- list(
  make_option("--input", type = "character"),
  make_option("--clusters", type = "character"),
  make_option("--groups", type = "character"),
  make_option("--output-dir", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)

# Minimal loader duplicated here to avoid a cross-script source() dance under
# Rscript (no reliable sys.frame / script-path handling). Keep in lock-step
# with load_and_cluster.R's load_input(); if you add a format there, add it
# here too.
load_input <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    obj <- readRDS(path)
    if (!inherits(obj, "Seurat")) stop("not a Seurat object")
    return(obj)
  }
  if (ext == "h5ad") {
    if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
      stop("SeuratDisk required for .h5ad")
    }
    h5seurat <- sub("\\.h5ad$", ".h5seurat", path, ignore.case = TRUE)
    if (!file.exists(h5seurat)) {
      SeuratDisk::Convert(path, dest = "h5seurat", overwrite = FALSE)
    }
    return(SeuratDisk::LoadH5Seurat(h5seurat))
  }
  stop("unsupported extension: .", ext)
}

message("[aggregate_rna] loading ", opt$input)
obj <- load_input(opt$input)

if (!"RNA" %in% Assays(obj)) {
  stop("no RNA assay found; cannot aggregate raw counts. ",
       "If this is an ATAC-only object, rerun the parent with --atac-only.")
}
# Seurat v5 replaced slot= with layer=; fall back gracefully for v4 objects.
counts <- tryCatch(
  GetAssayData(obj, assay = "RNA", layer = "counts"),
  error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
)
if (!inherits(counts, "dgCMatrix") && !inherits(counts, "Matrix")) {
  counts <- as(counts, "CsparseMatrix")
}

clusters_df <- read.csv(opt$clusters, stringsAsFactors = FALSE,
                         check.names = FALSE)
stopifnot("barcode" %in% colnames(clusters_df),
          "seurat_cluster" %in% colnames(clusters_df))

# Index barcodes.
bc_to_col <- setNames(seq_len(ncol(counts)), colnames(counts))
groups_plan <- fromJSON(opt$groups, simplifyDataFrame = FALSE)

aggregate_one <- function(group_spec) {
  cl <- group_spec$cluster
  metadata <- group_spec$metadata   # named list/dict: {col: val, ...} or NULL
  sel <- as.character(clusters_df$seurat_cluster) == as.character(cl)
  if (!is.null(metadata) && length(metadata) > 0) {
    for (col in names(metadata)) {
      if (!(col %in% colnames(clusters_df))) {
        warning("[aggregate_rna] group ", group_spec$name,
                ": metadata col '", col, "' not in clusters.csv; skipping.")
        return(NULL)
      }
      sel <- sel &
             !is.na(clusters_df[[col]]) &
             as.character(clusters_df[[col]]) == as.character(metadata[[col]])
    }
  }
  # For coembed objects (output of coembed-construct or any merged
  # RNA+ATAC object): the meta.data 'assay' column flags which modality
  # each cell came from. ATAC-origin cells in such objects carry IMPUTED
  # RNA counts (from TransferData), not measured RNA — summing those
  # would inflate the pseudobulk with model-derived values. Restrict to
  # RNA-origin cells when the column is present. For non-coembed inputs
  # (the column is absent), this is a no-op.
  if ("assay" %in% colnames(clusters_df)) {
    sel <- sel & as.character(clusters_df$assay) == "RNA"
  }
  cells <- clusters_df$barcode[sel]
  cells <- intersect(cells, names(bc_to_col))
  if (length(cells) == 0) {
    warning("[aggregate_rna] group ", group_spec$name,
            ": no cells match; skipping.")
    return(NULL)
  }
  idx <- bc_to_col[cells]
  # Sum across cells -> gene x 1 vector. rowSums on a sparse slice is the
  # fastest Matrix-pkg path; materializing a dense slice on 30k genes x 10k
  # cells is memory-expensive.
  summed <- Matrix::rowSums(counts[, idx, drop = FALSE])
  stopifnot(length(summed) == nrow(counts))
  out_path <- file.path(opt$`output-dir`, paste0(group_spec$name, ".tsv"))
  # GeneQuant format: gene_id<TAB>count, no header, no row names as a
  # separate column.
  write.table(
    data.frame(gene_id = rownames(counts), count = as.integer(summed)),
    file = out_path, sep = "\t",
    col.names = FALSE, row.names = FALSE, quote = FALSE
  )
  out_path
}

written <- list()
for (g in groups_plan$groups) {
  out <- aggregate_one(g)
  if (!is.null(out)) {
    written[[length(written) + 1]] <- list(name = g$name, path = out)
  }
}

# Write an index so the Python orchestrator and human reviewers can see what
# landed where. Not strictly necessary for the build-taiji-input handoff —
# the manifest.tsv already carries this info.
writeLines(
  toJSON(written, pretty = TRUE, auto_unbox = TRUE),
  file.path(opt$`output-dir`, "_index.json")
)

message("[aggregate_rna] wrote ", length(written), " RNA pseudobulk TSVs.")
