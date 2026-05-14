#!/usr/bin/env Rscript
# chooseR_clustering.R
#
# Resolution selection for co-embedded scRNA+scATAC objects via silhouette
# analysis across a resolution sweep.
#
# Algorithm (adapts chooseR, Ren et al. 2021):
#   1. FindNeighbors on joint PCA (if SNN not already in object)
#   2. For each candidate resolution:
#      a. FindClusters (Louvain)
#      b. Subsample --n-subsample cells uniformly
#      c. Compute Euclidean distances in PCA embedding (--n-pcs dims)
#      d. cluster::silhouette() on the subsample
#      e. Record mean silhouette width, n_clusters, median cluster size
#   3. Detect leading spike (sil[1] > CHOOSER_SPIKE_FACTOR * sil[2]); exclude
#      from selection. Find longest stable plateau (|Δsil| < CHOOSER_PLATEAU_TOL).
#      Select MIDPOINT of the plateau (not onset) — ensures the selected
#      resolution is clearly inside the stable zone, equidistant from both edges.
#      Fallback: highest silhouette among non-spike candidates.
#   4. Re-cluster at chosen resolution, drop clusters < --min-cluster-cells
#   5. Save updated object + chooseR_summary.json + silhouette plot + UMAP

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(cluster)
  library(jsonlite)
})

# -------------------------------------------------------------------------
# CLI
# -------------------------------------------------------------------------
option_list <- list(
  make_option("--input",            type = "character",
              help = "Path to coembed .rds (output of coembed.R with --no-cluster)"),
  make_option("--output",           type = "character",
              help = "Path to write the updated clustered .rds"),
  make_option("--resolutions",      type = "character",
              default = "0.1,0.2,0.4,0.6,0.8,1.0,1.2,1.5,2.0,2.5,3.0,4.0,5.0",
              help = "Comma-separated list of resolutions to sweep"),
  make_option("--n-subsample",      type = "integer", default = 2000L,
              help = "Cells sampled per resolution for silhouette (full 81k infeasible)"),
  make_option("--n-pcs",            type = "integer", default = 30L,
              help = "PCA dims to use for distance computation and FindNeighbors"),
  make_option("--min-cluster-cells", type = "integer", default = 20L,
              help = "Drop clusters below this size after final clustering"),
  make_option("--pca-reduction",    type = "character", default = "pca",
              help = "Name of the joint PCA reduction in the object"),
  make_option("--metadata-cols",    type = "character", default = NULL,
              help = "Comma-separated metadata cols to carry into UMAP panels"),
  make_option("--n-bootstraps",     type = "integer", default = 10L,
              help = "Bootstrap draws per resolution for silhouette stability estimate. Each draw independently resamples --n-subsample cells; mean±SD across draws is reported. Use 1 to disable bootstrapping (single draw, original behaviour)."),
  make_option("--seed",             type = "integer", default = 42L),
  make_option("--no-plot",          action = "store_true", default = FALSE)
)
opt <- parse_args(OptionParser(option_list = option_list))
for (k in c("input", "output")) {
  if (is.null(opt[[k]])) stop("missing required --", k)
}

set.seed(opt$seed)
output_dir <- dirname(opt$output)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

resolutions <- as.numeric(trimws(strsplit(opt$resolutions, ",", fixed = TRUE)[[1]]))
resolutions <- sort(unique(resolutions[!is.na(resolutions) & resolutions > 0]))
meta_cols <- if (!is.null(opt$`metadata-cols`)) {
  trimws(strsplit(opt$`metadata-cols`, ",", fixed = TRUE)[[1]])
} else character(0)

message("[chooseR] loading coembed from ", opt$input)
obj <- readRDS(opt$input)
message(sprintf("[chooseR]   %d cells, reductions: {%s}",
                ncol(obj), paste(names(obj@reductions), collapse = ",")))

# -------------------------------------------------------------------------
# FindNeighbors if SNN graph not present
# -------------------------------------------------------------------------
snn_graphs <- grep("_snn$", names(obj@graphs), value = TRUE)
if (length(snn_graphs) == 0) {
  message(sprintf("[chooseR] FindNeighbors on '%s' dims 1:%d",
                  opt$`pca-reduction`, opt$`n-pcs`))
  obj <- FindNeighbors(obj, reduction = opt$`pca-reduction`,
                       dims = 1:opt$`n-pcs`, verbose = FALSE)
  snn_graphs <- grep("_snn$", names(obj@graphs), value = TRUE)
  message("[chooseR]   SNN graph: ", paste(snn_graphs, collapse = ", "))
} else {
  message("[chooseR]   existing SNN graph found: ", paste(snn_graphs, collapse = ", "),
          " — skipping FindNeighbors")
}

# -------------------------------------------------------------------------
# Extract PCA embedding for silhouette distances
# -------------------------------------------------------------------------
if (!opt$`pca-reduction` %in% names(obj@reductions)) {
  stop("'", opt$`pca-reduction`, "' not found in object. Available: ",
       paste(names(obj@reductions), collapse = ", "))
}
pca_emb <- Embeddings(obj, reduction = opt$`pca-reduction`)[, 1:opt$`n-pcs`]
n_cells  <- nrow(pca_emb)
message(sprintf("[chooseR] PCA embedding: %d cells x %d dims", n_cells, opt$`n-pcs`))

# -------------------------------------------------------------------------
# Silhouette sweep
# -------------------------------------------------------------------------
message(sprintf("[chooseR] sweeping %d resolutions: %s",
                length(resolutions), paste(sprintf("%.1f", resolutions), collapse = ", ")))

results <- vector("list", length(resolutions))

for (i in seq_along(resolutions)) {
  r <- resolutions[i]
  obj_r <- FindClusters(obj, resolution = r, verbose = FALSE)
  clust  <- as.integer(Idents(obj_r))
  n_clust <- length(unique(clust))
  sizes  <- table(clust)
  med_size <- median(as.integer(sizes))

  if (n_clust < 2) {
    message(sprintf("[chooseR]   res=%.2f -> %d cluster(s) — skipping silhouette (need >=2)",
                    r, n_clust))
    results[[i]] <- list(resolution = r, n_clusters = n_clust,
                         median_cluster_size = as.numeric(med_size),
                         mean_silhouette = NA_real_,
                         n_subsample = 0L)
    next
  }

  # Bootstrap silhouette: repeat --n-bootstraps independent subsamples of
  # --n-subsample cells. FindClusters runs once per resolution (expensive);
  # the subsample + silhouette step is cheap and runs B times. Mean±SD across
  # draws is the stability proxy — wide SD relative to inter-resolution
  # differences means the subsample size is too small.
  n_boots   <- opt$`n-bootstraps`
  boot_sils <- rep(NA_real_, n_boots)

  for (b in seq_len(n_boots)) {
    idx <- sample(n_cells, min(opt$`n-subsample`, n_cells))
    sub_emb   <- pca_emb[idx, ]
    sub_clust <- clust[idx]

    sub_tab     <- table(sub_clust)
    multi_clust <- as.integer(names(sub_tab)[sub_tab >= 2])
    valid_idx   <- which(sub_clust %in% multi_clust)

    if (length(valid_idx) < 10 || length(unique(sub_clust[valid_idx])) < 2) next

    sub_emb_v   <- sub_emb[valid_idx, ]
    sub_clust_v <- sub_clust[valid_idx]

    d     <- dist(sub_emb_v)
    sil_b <- tryCatch(cluster::silhouette(sub_clust_v, d),
                      error = function(e) NULL)
    if (!is.null(sil_b)) boot_sils[b] <- mean(sil_b[, 3])
  }

  valid_boots <- boot_sils[!is.na(boot_sils)]
  if (length(valid_boots) == 0) {
    message(sprintf("[chooseR]   res=%.2f -> %d clusters, all bootstrap draws too thin — skipping",
                    r, n_clust))
    results[[i]] <- list(resolution = r, n_clusters = n_clust,
                         median_cluster_size = as.numeric(med_size),
                         mean_silhouette = NA_real_, sd_silhouette = NA_real_,
                         n_bootstraps = 0L, n_subsample = opt$`n-subsample`)
    next
  }

  mean_sil <- mean(valid_boots)
  sd_sil   <- if (length(valid_boots) > 1) sd(valid_boots) else 0

  message(sprintf("[chooseR]   res=%.2f -> %d clusters, median_size=%.0f, mean_sil=%.4f ± %.4f (n_boots=%d)",
                  r, n_clust, med_size, mean_sil, sd_sil, length(valid_boots)))
  results[[i]] <- list(resolution = r, n_clusters = n_clust,
                       median_cluster_size = as.numeric(med_size),
                       mean_silhouette = mean_sil, sd_silhouette = sd_sil,
                       n_bootstraps = length(valid_boots),
                       n_subsample = opt$`n-subsample`)
}

# -------------------------------------------------------------------------
# Select optimal resolution — plateau onset, not global maximum
# -------------------------------------------------------------------------
# The silhouette curve typically has a large spike at very low resolutions
# (where Louvain just separates the dominant broad compartments) and then
# a stable plateau at intermediate resolutions before collapsing to ~0 at
# high resolutions. The biologically meaningful choice is the ONSET of the
# plateau — the lowest resolution where the curve has stabilised — not the
# global maximum, which simply identifies the coarsest split.
#
# Algorithm:
#   1. Exclude the global maximum if it is an isolated spike: defined as
#      sil[1] > spike_factor * sil[2] (first point much higher than second).
#      When true, the spike is flagged and selection proceeds on the
#      remaining points.
#   2. Among the remaining valid points, find the plateau: the run of
#      consecutive resolutions where |sil[i] - sil[i-1]| < plateau_tol.
#      Select the FIRST (lowest resolution) point in the longest such run.
#   3. If no plateau is found (monotonically decreasing), fall back to the
#      valid point with the highest silhouette excluding the spike.
sil_vals   <- sapply(results, function(x) x$mean_silhouette)
valid_mask <- !is.na(sil_vals)
if (!any(valid_mask)) {
  stop("[chooseR] no valid silhouette scores computed — check --n-pcs and object layout")
}

spike_factor <- as.numeric(Sys.getenv("CHOOSER_SPIKE_FACTOR", unset = "1.5"))
plateau_tol  <- as.numeric(Sys.getenv("CHOOSER_PLATEAU_TOL",  unset = "0.05"))

valid_idx  <- which(valid_mask)
valid_sils <- sil_vals[valid_idx]

# Detect leading spike
spike_flagged <- FALSE
search_idx <- valid_idx
if (length(valid_sils) >= 2 &&
    valid_sils[1] > spike_factor * valid_sils[2]) {
  message(sprintf(
    "[chooseR] leading spike detected: res=%.2f sil=%.4f is >%.1fx res=%.2f sil=%.4f — excluding from plateau search",
    results[[valid_idx[1]]]$resolution, valid_sils[1],
    spike_factor,
    results[[valid_idx[2]]]$resolution, valid_sils[2]))
  spike_flagged <- TRUE
  search_idx  <- valid_idx[-1]
}

plateau_onset_idx <- NA_integer_
if (length(search_idx) >= 2) {
  search_sils <- sil_vals[search_idx]
  diffs <- abs(diff(search_sils))
  # Find runs of consecutive small diffs (plateau)
  in_plateau <- c(FALSE, diffs < plateau_tol)   # first point always starts a run
  best_run_len <- 0L
  best_run_start <- NA_integer_
  run_len <- 0L; run_start <- NA_integer_
  for (j in seq_along(search_idx)) {
    if (in_plateau[j]) {
      if (run_len == 0L) run_start <- j - 1L   # include the point before the small diff
      run_len <- run_len + 1L
    } else {
      if (run_len > best_run_len) {
        best_run_len   <- run_len
        best_run_start <- run_start
      }
      run_len <- 0L; run_start <- NA_integer_
    }
  }
  if (run_len > best_run_len) {
    best_run_len <- run_len; best_run_start <- run_start
  }
  if (!is.na(best_run_start) && best_run_len >= 1L) {
    # Select midpoint of the plateau run: equidistant from the leading edge
    # (where the curve may still be descending into the plateau) and the
    # trailing edge (where collapse begins). First point of run = best_run_start.
    plateau_start_pos <- max(1L, best_run_start)
    plateau_mid_pos   <- plateau_start_pos + floor(best_run_len / 2)
    plateau_onset_idx <- search_idx[min(plateau_mid_pos, length(search_idx))]
    message(sprintf(
      "[chooseR] plateau detected: run of %d resolutions (res=%.2f to res=%.2f); midpoint selected at res=%.2f (sil=%.4f)",
      best_run_len,
      results[[search_idx[plateau_start_pos]]]$resolution,
      results[[search_idx[min(plateau_start_pos + best_run_len, length(search_idx))]]]$resolution,
      results[[plateau_onset_idx]]$resolution,
      sil_vals[plateau_onset_idx]))
  }
}

if (!is.na(plateau_onset_idx)) {
  best_idx <- plateau_onset_idx
} else {
  # Fallback: highest silhouette among non-spike candidates
  best_idx <- search_idx[which.max(sil_vals[search_idx])]
  message("[chooseR] no stable plateau found; falling back to highest sil among non-spike candidates")
}

best_res <- results[[best_idx]]$resolution
best_sil <- sil_vals[best_idx]
message(sprintf("[chooseR] selected resolution: %.2f (mean_sil=%.4f, n_clusters=%d, spike_flagged=%s)",
                best_res, best_sil, results[[best_idx]]$n_clusters,
                spike_flagged))

# -------------------------------------------------------------------------
# Final clustering at chosen resolution
# -------------------------------------------------------------------------
message("[chooseR] final FindClusters at resolution=", best_res)
obj <- FindClusters(obj, resolution = best_res, verbose = FALSE)

cluster_vec <- as.character(Idents(obj))
sizes <- table(cluster_vec)
keep  <- names(sizes)[sizes >= opt$`min-cluster-cells`]
dropped <- setdiff(names(sizes), keep)
if (length(dropped) > 0) {
  message("[chooseR] dropping ", length(dropped), " clusters below --min-cluster-cells=",
          opt$`min-cluster-cells`, ": ", paste(dropped, collapse = ", "))
  obj <- subset(obj, cells = colnames(obj)[cluster_vec %in% keep])
}

n_kept_cells <- ncol(obj)
message(sprintf("[chooseR] final: %d clusters, %d cells (dropped %d small clusters)",
                length(keep), n_kept_cells, length(dropped)))

# -------------------------------------------------------------------------
# Per-cluster assay composition (RNA vs ATAC cells per cluster)
# -------------------------------------------------------------------------
if ("assay" %in% colnames(obj@meta.data)) {
  tab <- as.data.frame(table(obj$seurat_clusters, obj$assay))
  colnames(tab) <- c("cluster", "assay", "n_cells")
  # Warn on clusters with very few ATAC cells (will produce shallow peaks)
  atac_per_clust <- tab[tab$assay == "ATAC", ]
  thin_atac <- atac_per_clust[atac_per_clust$n_cells < 10, ]
  if (nrow(thin_atac) > 0) {
    message(sprintf("[chooseR] WARN: %d clusters have <10 ATAC cells — peak calls will be noisy: %s",
                    nrow(thin_atac), paste(thin_atac$cluster, collapse = ", ")))
  }
} else {
  tab <- NULL
}

# -------------------------------------------------------------------------
# Save object
# -------------------------------------------------------------------------
saveRDS(obj, opt$output)
message("[chooseR] saved clustered object -> ", opt$output)

# -------------------------------------------------------------------------
# Summary JSON
# -------------------------------------------------------------------------
summary_out <- list(
  input      = normalizePath(opt$input, mustWork = FALSE),
  output     = normalizePath(opt$output, mustWork = FALSE),
  n_cells    = n_kept_cells,
  n_pcs      = opt$`n-pcs`,
  n_subsample = opt$`n-subsample`,
  seed       = opt$seed,
  resolutions_swept = resolutions,
  sweep_results = results,
  chosen_resolution = best_res,
  chosen_mean_silhouette = best_sil,
  chosen_sd_silhouette = if (is.null(results[[best_idx]]$sd_silhouette)) NA else results[[best_idx]]$sd_silhouette,
  n_bootstraps  = opt$`n-bootstraps`,
  spike_flagged = spike_flagged,
  spike_factor  = spike_factor,
  plateau_tol   = plateau_tol,
  n_clusters_final = length(keep),
  dropped_clusters = dropped,
  min_cluster_cells = opt$`min-cluster-cells`
)
writeLines(
  toJSON(summary_out, pretty = TRUE, auto_unbox = TRUE, na = "null"),
  file.path(output_dir, "chooseR_summary.json")
)
message("[chooseR] summary -> ", file.path(output_dir, "chooseR_summary.json"))

# -------------------------------------------------------------------------
# Plots
# -------------------------------------------------------------------------
if (!opt$`no-plot` && requireNamespace("ggplot2", quietly = TRUE)) {
  qc_dir <- file.path(output_dir, "qc")
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

  # 1. Silhouette vs resolution
  sil_df <- data.frame(
    resolution   = sapply(results, `[[`, "resolution"),
    mean_sil     = sapply(results, `[[`, "mean_silhouette"),
    sd_sil       = sapply(results, function(x) if (is.null(x$sd_silhouette)) NA_real_ else x$sd_silhouette),
    n_clusters   = sapply(results, `[[`, "n_clusters"),
    stringsAsFactors = FALSE
  )
  sil_df <- sil_df[!is.na(sil_df$mean_sil), ]

  # Mark the spike (excluded from selection) in grey if it was flagged
  spike_res <- if (spike_flagged) results[[valid_idx[1]]]$resolution else NA_real_
  p_sil <- ggplot2::ggplot(sil_df,
    ggplot2::aes(x = resolution, y = mean_sil)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_sil - sd_sil, ymax = mean_sil + sd_sil),
      width = 0.04, colour = "#aec7e8", linewidth = 0.6, na.rm = TRUE) +
    ggplot2::geom_line(colour = "#1f77b4", linewidth = 0.8) +
    ggplot2::geom_point(size = 2.5, colour = "#1f77b4")
  if (!is.na(spike_res)) {
    p_sil <- p_sil +
      ggplot2::geom_point(data = sil_df[sil_df$resolution == spike_res, ],
                          ggplot2::aes(x = resolution, y = mean_sil),
                          colour = "#aaaaaa", size = 3.5, shape = 4) +
      ggplot2::annotate("text", x = spike_res,
                        y = sil_df$mean_sil[sil_df$resolution == spike_res],
                        label = "spike\n(excluded)", hjust = -0.1,
                        colour = "#aaaaaa", size = 2.8)
  }
  p_sil <- p_sil +
    ggplot2::geom_vline(xintercept = best_res, linetype = "dashed",
                        colour = "#d62728", linewidth = 0.7) +
    ggplot2::annotate("text", x = best_res, y = min(sil_df$mean_sil),
                      label = sprintf("plateau midpoint\nres=%.2f", best_res),
                      hjust = -0.1, colour = "#d62728", size = 3) +
    ggplot2::labs(
      title    = sprintf("chooseR: mean silhouette vs. resolution (n_sub=%d)", opt$`n-subsample`),
      subtitle = sprintf("plateau midpoint res=%.2f (sil=%.3f±%.3f, %d clusters, %d bootstraps)%s",
                         best_res, best_sil,
                         if (is.null(results[[best_idx]]$sd_silhouette) || is.na(results[[best_idx]]$sd_silhouette)) 0 else results[[best_idx]]$sd_silhouette,
                         results[[best_idx]]$n_clusters, opt$`n-bootstraps`,
                         if (spike_flagged) " — leading spike excluded" else ""),
      x = "Louvain resolution", y = "mean silhouette width") +
    ggplot2::theme_bw(base_size = 12)

  ggplot2::ggsave(file.path(qc_dir, "chooseR_silhouette.png"),
                  plot = p_sil, width = 7, height = 4, dpi = 120, bg = "white")
  message("[chooseR] silhouette plot -> ", file.path(qc_dir, "chooseR_silhouette.png"))

  # 2. n_clusters vs resolution (secondary diagnostic)
  sil_df_all <- data.frame(
    resolution = sapply(results, `[[`, "resolution"),
    n_clusters = sapply(results, `[[`, "n_clusters"),
    stringsAsFactors = FALSE
  )
  p_nclust <- ggplot2::ggplot(sil_df_all,
    ggplot2::aes(x = resolution, y = n_clusters)) +
    ggplot2::geom_line(colour = "#2ca02c", linewidth = 0.8) +
    ggplot2::geom_point(size = 2.5, colour = "#2ca02c") +
    ggplot2::geom_vline(xintercept = best_res, linetype = "dashed",
                        colour = "#d62728", linewidth = 0.7) +
    ggplot2::labs(title = "Number of clusters vs. resolution",
                  x = "Louvain resolution", y = "n clusters") +
    ggplot2::theme_bw(base_size = 12)

  ggplot2::ggsave(file.path(qc_dir, "chooseR_nclusters.png"),
                  plot = p_nclust, width = 7, height = 4, dpi = 120, bg = "white")

  # 3. UMAP colored by final clusters + metadata cols
  if ("umap" %in% names(obj@reductions)) {
    groupby_cols <- c("seurat_clusters", "assay",
                      intersect(meta_cols, colnames(obj@meta.data)))
    panels <- list()
    for (col in groupby_cols) {
      p <- tryCatch(
        DimPlot(obj, reduction = "umap", group.by = col,
                label = (col == "seurat_clusters"), repel = TRUE) +
          ggplot2::ggtitle(col),
        error = function(e) NULL
      )
      if (!is.null(p)) panels[[length(panels) + 1]] <- p
    }
    if (length(panels) > 0 && requireNamespace("patchwork", quietly = TRUE)) {
      ncol_p <- min(2, length(panels))
      n_rows  <- ceiling(length(panels) / ncol_p)
      ggplot2::ggsave(
        file.path(qc_dir, "umap_chooseR.png"),
        plot = patchwork::wrap_plots(panels, ncol = ncol_p),
        width = 6 * ncol_p, height = 5 * n_rows, dpi = 100, bg = "white"
      )
      message("[chooseR] UMAP -> ", file.path(qc_dir, "umap_chooseR.png"))
    }
  }
}

message("[chooseR] done.")
