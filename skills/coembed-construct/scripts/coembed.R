#!/usr/bin/env Rscript
# coembed.R
#
# Co-embed scRNA-seq + scATAC-seq (separate-assay) into a shared latent
# space. Implements the Stuart 2019 / Signac integrate_atac vignette
# end-to-end:
#
#   1. RNA: NormalizeData -> FindVariableFeatures -> ScaleData -> RunPCA -> RunUMAP
#   2. ATAC: RunTFIDF -> FindTopFeatures -> RunSVD -> RunUMAP
#   3. Set EnsDb annotation on ATAC -> GeneActivity -> ACTIVITY assay
#   4. FindTransferAnchors (RNA reference, ATAC query, CCA reduction)
#   5. TransferData on RNA expression matrix -> imputed RNA assay on ATAC
#   6. merge(rna, atac) -> coembed
#   7. ScaleData + RunPCA + RunUMAP on the merged object (RNA features)
#   8. FindNeighbors + FindClusters with scale-aware resolution binary search
#   9. saveRDS + qc/umap.png + coembed_summary.json
#
# Output coembed.rds is drop-in for pseudobulk-construct's
# --use-existing-clusters mode. Cells are tagged via meta.data$assay
# ("RNA" or "ATAC") so downstream skills can filter origin-aware:
# aggregate_rna only sums real (assay=='RNA') counts, call_peaks only
# scopes fragments to ATAC-origin cells.
#
# The skill is deliberately narrow: it does NOT do cell-type label
# transfer (TransferData(refdata = seurat_annotations) from the vignette).
# Cluster IDs come from de novo clustering on the shared space; named
# labels are out of scope (see SKILL.md for rationale).

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(jsonlite)
})

# ------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------

option_list <- list(
  make_option("--rna",                type = "character",
              help = "Path to scRNA-seq Seurat .rds"),
  make_option("--atac",               type = "character",
              help = "Path to scATAC-seq Seurat .rds (ChromatinAssay required)"),
  make_option("--output",             type = "character",
              help = "Path to write coembed.rds"),
  make_option("--genome",             type = "character",
              help = "hg38 / hg19 / mm10 / mm39"),
  make_option("--target-cluster-size", type = "integer", default = 200L),
  make_option("--min-cluster-cells",   type = "integer", default = 20L),
  make_option("--metadata-cols",       type = "character", default = NULL,
              help = "Comma-separated metadata cols to preserve in the merged object."),
  make_option("--n-pcs",               type = "integer", default = 30L,
              help = "Joint PCA dims for FindNeighbors / RunUMAP."),
  make_option("--lsi-skip-first",      type = "integer", default = 1L,
              help = "Skip first N LSI components on ATAC (depth-correlated)."),
  make_option("--resolution",          type = "double", default = NA_real_,
              help = "Force a single resolution; default: scale-aware binary search."),
  make_option("--no-cluster",          action = "store_true", default = FALSE),
  make_option("--no-plot",             action = "store_true", default = FALSE)
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("rna", "atac", "output", "genome")
for (k in required) if (is.null(opt[[k]])) stop("missing required --", k)

output_dir <- dirname(opt$output)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

meta_cols <- if (!is.null(opt$`metadata-cols`)) {
  trimws(strsplit(opt$`metadata-cols`, ",", fixed = TRUE)[[1]])
} else character(0)

# ------------------------------------------------------------------------
# EnsDb dispatch (for GeneActivity annotation)
# ------------------------------------------------------------------------

load_ensdb <- function(genome) {
  pkg <- switch(tolower(genome),
                hg38 = "EnsDb.Hsapiens.v86",
                hg19 = "EnsDb.Hsapiens.v86",
                mm10 = "EnsDb.Mmusculus.v79",
                mm39 = "EnsDb.Mmusculus.v79",
                NULL)
  if (is.null(pkg)) {
    stop("unsupported --genome '", genome,
         "'. Add an EnsDb dispatch for it in coembed.R.")
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(pkg, " is required for genome ", genome,
         ". Install via Bioconductor: BiocManager::install('", pkg, "').")
  }
  getFromNamespace(pkg, pkg)
}

# ------------------------------------------------------------------------
# Pick the ATAC assay from a Seurat object
# ------------------------------------------------------------------------

pick_atac_assay <- function(obj) {
  for (name in c("ATAC", "peaks")) {
    if (name %in% Assays(obj)) return(name)
  }
  hits <- Assays(obj)[sapply(Assays(obj),
                             function(a) inherits(obj[[a]], "ChromatinAssay"))]
  if (length(hits) == 0) {
    stop("no ATAC/ChromatinAssay found in ATAC object. ",
         "Available assays: ", paste(Assays(obj), collapse = ", "))
  }
  hits[1]
}

# ------------------------------------------------------------------------
# Resolution binary search (adapted from load_and_cluster.R; targets a
# given mean cluster size on the JOINT FindNeighbors graph)
# ------------------------------------------------------------------------

seed_resolution <- function(n_cells, target_size) {
  target_n_clusters <- max(2, n_cells / target_size)
  r0 <- 0.15 * log2(max(2, n_cells / 1000))
  min(max(r0 * sqrt(target_n_clusters / 20), 0.05), 3.0)
}

in_range <- function(mean_size, target_size) {
  if (isTRUE(all.equal(target_size, 200))) {
    return(mean_size >= 100 && mean_size <= 300)
  }
  lo <- 0.5 * target_size
  hi <- 1.5 * target_size
  mean_size >= lo && mean_size <= hi
}

binary_search_resolution <- function(obj, n_cells, target_size, max_iter = 8L) {
  r <- seed_resolution(n_cells, target_size)
  lo <- 0.01; hi <- 5.0
  trace <- list()
  obj_try <- obj
  for (i in seq_len(max_iter)) {
    obj_try <- FindClusters(obj, resolution = r, verbose = FALSE)
    n_clust <- length(unique(Idents(obj_try)))
    mean_size <- n_cells / n_clust
    trace[[length(trace) + 1]] <- list(
      iter = i, resolution = r,
      n_clusters = as.integer(n_clust),
      mean_cluster_size = round(mean_size, 1)
    )
    message(sprintf(
      "[coembed] iter %d: res=%.3f -> %d clusters, mean_size=%.1f",
      i, r, n_clust, mean_size))
    if (in_range(mean_size, target_size)) {
      return(list(obj = obj_try, resolution = r, trace = trace))
    }
    if (mean_size > target_size) {
      lo <- r
      r <- if (hi < 5.0) 0.5 * (r + hi) else min(r * 1.5, hi)
    } else {
      hi <- r
      r <- if (lo > 0.01) 0.5 * (r + lo) else max(r * 0.67, lo)
    }
  }
  warning("[coembed] resolution search did not converge; ",
          "returning last iteration.")
  list(obj = obj_try, resolution = r, trace = trace)
}

# ------------------------------------------------------------------------
# Stage 1 — load & process RNA
# ------------------------------------------------------------------------

message("[coembed] loading RNA from ", opt$rna)
rna <- readRDS(opt$rna)
if (!inherits(rna, "Seurat")) stop("RNA file is not a Seurat object")
DefaultAssay(rna) <- "RNA"

if (is.null(GetAssayData(rna, slot = "counts")) ||
    nrow(GetAssayData(rna, slot = "counts")) == 0) {
  stop("RNA object has no raw counts. Cannot proceed with NormalizeData.")
}

message("[coembed] RNA: standard preprocessing (Norm/VF/Scale/PCA/UMAP)")
rna <- NormalizeData(rna, verbose = FALSE)
rna <- FindVariableFeatures(rna, verbose = FALSE)
rna <- ScaleData(rna, verbose = FALSE)
rna <- RunPCA(rna, npcs = opt$`n-pcs`, verbose = FALSE)
rna <- RunUMAP(rna, dims = 1:opt$`n-pcs`, verbose = FALSE)
rna$assay <- "RNA"

# ------------------------------------------------------------------------
# Stage 2 — load & process ATAC, set annotation, RunTFIDF/SVD/UMAP
# ------------------------------------------------------------------------

message("[coembed] loading ATAC from ", opt$atac)
atac <- readRDS(opt$atac)
if (!inherits(atac, "Seurat")) stop("ATAC file is not a Seurat object")

atac_assay <- pick_atac_assay(atac)
DefaultAssay(atac) <- atac_assay
message("[coembed] ATAC: using assay '", atac_assay, "'")

# Annotation pull for GeneActivity. Signac wants UCSC-style chromosome names
# and a `genome` attribute on the GRanges; this is what the vignette does.
ensdb <- load_ensdb(opt$genome)
message("[coembed] pulling gene annotations from EnsDb (", opt$genome, ") ...")
annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- opt$genome
Annotation(atac) <- annotations

message("[coembed] ATAC: TF-IDF / FindTopFeatures / RunSVD / RunUMAP")
atac <- RunTFIDF(atac, verbose = FALSE)
atac <- FindTopFeatures(atac, min.cutoff = "q0", verbose = FALSE)
atac <- RunSVD(atac, verbose = FALSE)
lsi_dims <- (opt$`lsi-skip-first` + 1):opt$`n-pcs`
atac <- RunUMAP(atac, reduction = "lsi", dims = lsi_dims,
                reduction.name = "umap.atac", verbose = FALSE)
atac$assay <- "ATAC"

# ------------------------------------------------------------------------
# Stage 3 — GeneActivity -> ACTIVITY assay
# ------------------------------------------------------------------------

message("[coembed] computing per-cell gene activities (ATAC -> ACTIVITY assay)")
gene_activities <- GeneActivity(atac, features = VariableFeatures(rna))
atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene_activities)
DefaultAssay(atac) <- "ACTIVITY"
atac <- NormalizeData(atac, verbose = FALSE)
atac <- ScaleData(atac, features = rownames(atac), verbose = FALSE)

# ------------------------------------------------------------------------
# Stage 4 — FindTransferAnchors (RNA ref, ATAC query, CCA)
# ------------------------------------------------------------------------

message("[coembed] FindTransferAnchors: RNA (ref) -> ATAC (query, ACTIVITY assay), CCA")
transfer_anchors <- FindTransferAnchors(
  reference       = rna,
  query           = atac,
  features        = VariableFeatures(rna),
  reference.assay = "RNA",
  query.assay     = "ACTIVITY",
  reduction       = "cca",
  verbose         = FALSE
)
n_anchors <- nrow(slot(transfer_anchors, "anchors"))
message("[coembed]   found ", n_anchors, " anchors")

# ------------------------------------------------------------------------
# Stage 5 — Impute RNA into ATAC cells via TransferData
# ------------------------------------------------------------------------

message("[coembed] TransferData: imputing RNA expression into ATAC cells")
genes_use <- VariableFeatures(rna)
refdata <- GetAssayData(rna, assay = "RNA", slot = "data")[genes_use, ]
imputation <- TransferData(
  anchorset        = transfer_anchors,
  refdata          = refdata,
  weight.reduction = atac[["lsi"]],
  dims             = lsi_dims,
  verbose          = FALSE
)
atac[["RNA"]] <- imputation

# ------------------------------------------------------------------------
# Stage 6 — Merge RNA + ATAC
# ------------------------------------------------------------------------

message("[coembed] merging RNA + ATAC objects")
DefaultAssay(rna)  <- "RNA"
DefaultAssay(atac) <- "RNA"
coembed <- merge(x = rna, y = atac)
message("[coembed]   merged: ", ncol(coembed), " cells (",
        sum(coembed$assay == "RNA"), " RNA + ",
        sum(coembed$assay == "ATAC"), " ATAC)")

# ------------------------------------------------------------------------
# Stage 7 — Joint PCA + UMAP on shared space (RNA features)
# ------------------------------------------------------------------------

message("[coembed] joint ScaleData / RunPCA / RunUMAP on shared space")
DefaultAssay(coembed) <- "RNA"
coembed <- ScaleData(coembed, features = genes_use, do.scale = FALSE, verbose = FALSE)
coembed <- RunPCA(coembed, features = genes_use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:opt$`n-pcs`, verbose = FALSE)

# ------------------------------------------------------------------------
# Stage 8 — Cluster on the shared space
# ------------------------------------------------------------------------

resolution_trace <- list()
chosen_resolution <- NA_real_
n_clusters_total <- NA_integer_
n_clusters_kept  <- NA_integer_
dropped_clusters <- character(0)

if (!opt$`no-cluster`) {
  message("[coembed] FindNeighbors on joint PCA")
  coembed <- FindNeighbors(coembed, dims = 1:opt$`n-pcs`, verbose = FALSE)

  if (!is.na(opt$resolution)) {
    message("[coembed] FindClusters at fixed resolution=", opt$resolution)
    coembed <- FindClusters(coembed, resolution = opt$resolution, verbose = FALSE)
    chosen_resolution <- opt$resolution
    resolution_trace <- list(list(
      iter = 1, resolution = opt$resolution,
      n_clusters = length(unique(Idents(coembed))),
      mean_cluster_size = ncol(coembed) / length(unique(Idents(coembed)))
    ))
  } else {
    res <- binary_search_resolution(coembed, ncol(coembed),
                                    opt$`target-cluster-size`)
    coembed <- res$obj
    chosen_resolution <- res$resolution
    resolution_trace <- res$trace
  }

  # Drop small clusters
  cluster_vec <- as.character(Idents(coembed))
  sizes <- table(cluster_vec)
  keep <- names(sizes)[sizes >= opt$`min-cluster-cells`]
  dropped_clusters <- setdiff(names(sizes), keep)
  if (length(dropped_clusters) > 0) {
    message("[coembed] dropping ", length(dropped_clusters),
            " clusters below --min-cluster-cells=", opt$`min-cluster-cells`,
            ": ", paste(dropped_clusters, collapse = ","))
    coembed <- subset(coembed, cells = colnames(coembed)[cluster_vec %in% keep])
  }
  n_clusters_total <- length(sizes)
  n_clusters_kept  <- length(keep)
} else {
  message("[coembed] --no-cluster set; skipping FindNeighbors + FindClusters")
}

# ------------------------------------------------------------------------
# Stage 9 — Save coembed object + summary JSON
# ------------------------------------------------------------------------

saveRDS(coembed, opt$output)
message("[coembed] wrote ", opt$output)

per_cluster_assay <- if (!opt$`no-cluster`) {
  tab <- as.data.frame(table(coembed$seurat_clusters, coembed$assay))
  colnames(tab) <- c("cluster", "assay", "n_cells")
  apply(tab, 1, as.list)
} else NULL

summary_json <- list(
  rna_input  = normalizePath(opt$rna, mustWork = FALSE),
  atac_input = normalizePath(opt$atac, mustWork = FALSE),
  genome     = opt$genome,
  n_cells_rna  = sum(coembed$assay == "RNA"),
  n_cells_atac = sum(coembed$assay == "ATAC"),
  n_cells_total = ncol(coembed),
  n_anchors  = n_anchors,
  n_genes_used = length(genes_use),
  chosen_resolution = chosen_resolution,
  resolution_trace = resolution_trace,
  n_clusters_total = n_clusters_total,
  n_clusters_kept  = n_clusters_kept,
  dropped_clusters = dropped_clusters,
  metadata_cols = meta_cols,
  output_path = normalizePath(opt$output, mustWork = FALSE)
)
writeLines(
  toJSON(summary_json, pretty = TRUE, auto_unbox = TRUE, na = "null"),
  file.path(output_dir, "coembed_summary.json")
)

# ------------------------------------------------------------------------
# QC UMAP — same panels load_and_cluster.R writes
# ------------------------------------------------------------------------

if (!opt$`no-plot` && requireNamespace("ggplot2", quietly = TRUE)) {
  qc_dir <- file.path(output_dir, "qc")
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

  emb <- Embeddings(coembed, reduction = "umap")
  umap_df <- data.frame(
    barcode = rownames(emb),
    umap_1  = emb[, 1],
    umap_2  = emb[, 2],
    assay   = coembed$assay,
    seurat_cluster = if (!opt$`no-cluster`) as.character(Idents(coembed)) else NA,
    stringsAsFactors = FALSE
  )
  write.csv(umap_df, file.path(qc_dir, "umap_coords.csv"), row.names = FALSE)

  groupby_cols <- c()
  if (!opt$`no-cluster`) groupby_cols <- c(groupby_cols, "seurat_clusters")
  groupby_cols <- c(groupby_cols, "assay")
  groupby_cols <- c(groupby_cols,
                    intersect(meta_cols, colnames(coembed@meta.data)))
  groupby_cols <- unique(groupby_cols)

  panels <- list()
  for (col in groupby_cols) {
    p <- tryCatch(
      DimPlot(coembed, reduction = "umap", group.by = col,
              label = (col == "seurat_clusters"), repel = TRUE) +
        ggplot2::ggtitle(col),
      error = function(e) {
        message("[coembed] DimPlot[", col, "] failed: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(p)) panels[[length(panels) + 1]] <- p
  }

  if (length(panels) > 0) {
    out_png <- file.path(qc_dir, "umap.png")
    if (requireNamespace("patchwork", quietly = TRUE) && length(panels) > 1) {
      ncol_plot <- min(2, length(panels))
      n_rows <- ceiling(length(panels) / ncol_plot)
      combined <- patchwork::wrap_plots(panels, ncol = ncol_plot)
      ggplot2::ggsave(out_png, plot = combined,
                      width = 6 * ncol_plot, height = 5 * n_rows,
                      dpi = 100, bg = "white")
    } else {
      ggplot2::ggsave(out_png, plot = panels[[1]],
                      width = 7, height = 5, dpi = 100, bg = "white")
    }
    message("[coembed] QC UMAP -> ", out_png)
  }
}

message("[coembed] done.")
