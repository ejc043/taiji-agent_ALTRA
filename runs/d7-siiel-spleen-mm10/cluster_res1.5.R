#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

res      <- 1.5
input    <- "runs/d7-siiel-spleen-mm10/coembed/coembed.rds"
out_rds  <- "runs/d7-siiel-spleen-mm10/coembed/coembed_clustered_res1.5.rds"
qc_dir   <- "runs/d7-siiel-spleen-mm10/coembed/qc"
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

message("loading ", input)
obj <- readRDS(input)
message(sprintf("  %d cells", ncol(obj)))

message("FindNeighbors (dims 1:30)")
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, verbose = FALSE)

message("FindClusters at res=", res)
obj <- FindClusters(obj, resolution = res, verbose = FALSE)
message(sprintf("  %d clusters", length(unique(Idents(obj)))))

saveRDS(obj, out_rds)
message("saved -> ", out_rds)

# UMAP panels: clusters, assay, tissue, genotype
groupby <- c("seurat_clusters", "assay", "tissue", "genotype")
groupby <- intersect(groupby, c("seurat_clusters", "assay", colnames(obj@meta.data)))

n_clust <- length(unique(Idents(obj)))
clust_cols <- rep_len(pals::alphabet(), n_clust)

panels <- lapply(groupby, function(col) {
  tryCatch({
    p <- DimPlot(obj, reduction = "umap", group.by = col,
                 label = (col == "seurat_clusters"), repel = TRUE,
                 cols = if (col == "seurat_clusters") clust_cols else NULL) +
      ggtitle(if (col == "seurat_clusters")
                sprintf("clusters (res=%.1f, n=%d)", res, n_clust)
              else col)
    p
  }, error = function(e) NULL)
})
panels <- Filter(Negate(is.null), panels)

out_png <- file.path(qc_dir, "umap_res1.5.png")
ncol_p  <- min(2, length(panels))
ggsave(out_png,
       plot = wrap_plots(panels, ncol = ncol_p),
       width = 6 * ncol_p, height = 5 * ceiling(length(panels) / ncol_p),
       dpi = 120, bg = "white")
message("UMAP -> ", out_png)
message("done.")
