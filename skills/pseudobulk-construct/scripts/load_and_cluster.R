#!/usr/bin/env Rscript
# load_and_cluster.R
#
# Load a single-cell object (.rds / .h5ad), coerce to Seurat, build a
# WNN / RNA-only / ATAC-only clustering with a scale-aware resolution search
# targeting ~target_cluster_size cells/cluster (range [100, 300] by default),
# drop clusters below the cell-count floor, auto-detect or accept metadata
# columns, and write:
#   - clusters.csv            (one row per cell: barcode, seurat_cluster, meta cols)
#   - resolution_trace.json   (every (resolution, n_clusters, mean_size) tried)
#   - per_cluster_barcodes/<cluster>.txt   (one barcode per line, for MACS2)
#   - groups_plan.json        (list of per-(cluster x metadata_col x value) triples)
#
# This is a single-shot CLI. It does no caching. The Python orchestrator
# handles dependency checks and gating; this script assumes Rscript, Seurat,
# Signac, optionally SeuratDisk are available.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(jsonlite)
  library(dplyr)
})

# ------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------

option_list <- list(
  make_option("--input", type = "character", help = "Path to .rds/.h5ad"),
  make_option("--output-dir", type = "character"),
  make_option("--signal", type = "character", default = "wnn",
              help = "wnn | rna | atac"),
  make_option("--target-cluster-size", type = "integer", default = 200L),
  make_option("--min-cluster-cells", type = "integer", default = 20L),
  make_option("--transferred-label-col", type = "character",
              default = "predicted.id"),
  make_option("--metadata-cols", type = "character", default = NULL,
              help = "Comma-separated columns to stratify on. Multiple columns produce the FULL CROSS-PRODUCT within each cluster (e.g. genotype,tissue -> per-cluster (WT,spleen) (WT,siiel) (KO,spleen) (KO,siiel))."),
  make_option("--cohort-col", type = "character", default = NULL,
              help = "Which metadata column's value becomes the manifest 'cohort' label. Default: first --metadata-cols entry."),
  make_option("--use-existing-clusters", action = "store_true", default = FALSE,
              help = "Skip clustering. Read seurat_clusters from the input object's @meta.data and use the existing UMAP reduction (or compute one from the existing PCA). Use this when the input is a coembed-construct output that's already been clustered on a shared embedding — re-clustering on a single modality's LSI/PCA would discard the joint structure."),
  make_option("--no-plot", action = "store_true", default = FALSE,
              help = "Skip the QC UMAP rendering step (default: write qc/umap.png)."),
  make_option("--yes", action = "store_true", default = FALSE)
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$input), !is.null(opt$`output-dir`))
output_dir <- opt$`output-dir`
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "per_cluster_barcodes"),
           recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------
# Object loading: dispatch on extension
# ------------------------------------------------------------------------

load_input <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    obj <- readRDS(path)
    if (!inherits(obj, "Seurat")) {
      stop("Loaded .rds is not a Seurat object (got ",
           paste(class(obj), collapse = "/"), "). Convert first.")
    }
    return(obj)
  }
  if (ext == "h5ad") {
    if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
      stop("SeuratDisk is required for .h5ad input. install.packages('SeuratDisk').")
    }
    # Convert to h5seurat side-by-side, then load. This is SeuratDisk's
    # canonical path; sceasy is an alternative but adds a Python dependency
    # we're trying to avoid at the R layer.
    h5seurat <- sub("\\.h5ad$", ".h5seurat", path, ignore.case = TRUE)
    if (!file.exists(h5seurat)) {
      SeuratDisk::Convert(path, dest = "h5seurat", overwrite = FALSE)
    }
    return(SeuratDisk::LoadH5Seurat(h5seurat))
  }
  stop("unsupported input extension: .", ext,
       " (expected .rds / .h5ad)")
}

message("[load_and_cluster] loading ", opt$input)
obj <- load_input(opt$input)
message("[load_and_cluster] loaded: ", ncol(obj), " cells, ",
        length(Assays(obj)), " assays (", paste(Assays(obj), collapse = ","), ")")

# ------------------------------------------------------------------------
# Assay preparation
# ------------------------------------------------------------------------

has_atac_assay <- function(obj) {
  any(c("ATAC", "peaks") %in% Assays(obj)) ||
    any(sapply(Assays(obj), function(a) inherits(obj[[a]], "ChromatinAssay")))
}

prep_rna <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
  obj
}

prep_atac <- function(obj) {
  atac_assay <- intersect(c("ATAC", "peaks"), Assays(obj))[1]
  if (is.na(atac_assay)) {
    # fall back: first ChromatinAssay
    hits <- Assays(obj)[sapply(Assays(obj),
                               function(a) inherits(obj[[a]], "ChromatinAssay"))]
    atac_assay <- hits[1]
  }
  if (is.na(atac_assay) || is.null(atac_assay)) {
    stop("no ATAC assay found; cannot cluster on ATAC LSI.")
  }
  DefaultAssay(obj) <- atac_assay
  # Signac canonical LSI pipeline.
  obj <- Signac::RunTFIDF(obj, verbose = FALSE)
  obj <- Signac::FindTopFeatures(obj, min.cutoff = "q5", verbose = FALSE)
  obj <- Signac::RunSVD(obj, verbose = FALSE)
  obj
}

# ------------------------------------------------------------------------
# Scale-aware resolution binary search
# ------------------------------------------------------------------------
#
# User guidance: larger datasets need higher resolution. We seed from
#   N_cells / target_cluster_size -> approximate number of clusters,
# then map target n_clusters to a starting Louvain resolution via a
# coarse heuristic (resolution ~ 0.15 * log2(N_cells / 1000), clamped to
# [0.05, 3.0]). Binary search from there, up to max_iter evaluations.
#
# The graph (SNN) is computed once; FindClusters just reruns the modularity
# step on each resolution, which is cheap.

seed_resolution <- function(n_cells, target_size) {
  target_n_clusters <- max(2, n_cells / target_size)
  r0 <- 0.15 * log2(max(2, n_cells / 1000))
  # Rough monotonic mapping: more target clusters -> higher resolution.
  # Clamp hard so we never start in pathological regimes.
  min(max(r0 * sqrt(target_n_clusters / 20), 0.05), 3.0)
}

in_range <- function(mean_size, target_size) {
  # The spec says ~100-300 cells/cluster. Center on target_size with a
  # symmetric 50% tolerance unless target is explicitly 200 (the spec
  # midpoint), in which case use the spec's own [100, 300] band.
  if (isTRUE(all.equal(target_size, 200))) {
    return(mean_size >= 100 && mean_size <= 300)
  }
  lo <- 0.5 * target_size
  hi <- 1.5 * target_size
  mean_size >= lo && mean_size <= hi
}

binary_search_resolution <- function(obj,
                                     graph_name,
                                     n_cells,
                                     target_size,
                                     max_iter = 8L) {
  r <- seed_resolution(n_cells, target_size)
  lo <- 0.01
  hi <- 5.0
  trace <- list()

  for (i in seq_len(max_iter)) {
    obj_try <- FindClusters(obj, graph.name = graph_name,
                            resolution = r, verbose = FALSE)
    clusters <- Idents(obj_try)
    n_clust <- length(unique(clusters))
    mean_size <- n_cells / n_clust
    trace[[length(trace) + 1]] <- list(
      iter = i, resolution = r,
      n_clusters = as.integer(n_clust),
      mean_cluster_size = round(mean_size, 1)
    )
    message(sprintf(
      "[load_and_cluster] iter %d: res=%.3f -> %d clusters, mean_size=%.1f",
      i, r, n_clust, mean_size))
    if (in_range(mean_size, target_size)) {
      return(list(obj = obj_try, resolution = r, trace = trace))
    }
    # Mean size too big -> need higher resolution; too small -> lower.
    if (mean_size > target_size) {
      lo <- r
      r <- if (hi < 5.0) 0.5 * (r + hi) else min(r * 1.5, hi)
    } else {
      hi <- r
      r <- if (lo > 0.01) 0.5 * (r + lo) else max(r * 0.67, lo)
    }
  }

  # Best-effort return: the last try is the closest-to-target we got.
  warning("[load_and_cluster] resolution search did not converge; ",
          "returning last iteration. Inspect resolution_trace.json.")
  list(obj = obj_try, resolution = r, trace = trace)
}

# ------------------------------------------------------------------------
# Run the chosen clustering signal
# ------------------------------------------------------------------------

n_cells <- ncol(obj)
message("[load_and_cluster] signal=", opt$signal, " on ", n_cells, " cells")

# ------------------------------------------------------------------------
# Branch A: Use clusters that already live in the object (coembed-construct
# output, or any user-pre-clustered .rds).
# ------------------------------------------------------------------------
if (opt$`use-existing-clusters`) {
  if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
    stop("--use-existing-clusters set but obj@meta.data has no ",
         "'seurat_clusters' column. Either cluster the object upstream, ",
         "or omit --use-existing-clusters to let this skill cluster.")
  }
  message("[load_and_cluster] using existing seurat_clusters from input ",
          "(--use-existing-clusters)")
  Idents(obj) <- "seurat_clusters"
  res_result <- list(
    obj = obj,
    resolution = NA_real_,
    trace = list(list(note = "skipped — used existing seurat_clusters"))
  )
  graph_name <- NA_character_
} else if (opt$signal == "wnn") {
  if (!has_atac_assay(obj)) {
    stop("WNN requested but no ATAC assay found in the object. ",
         "Rerun with --clustering-signal rna, or supply a multiome object.")
  }
  obj <- prep_rna(obj)
  obj <- prep_atac(obj)
  DefaultAssay(obj) <- "RNA"
  # Multi-modal neighbors. This function is the entire reason WNN exists:
  # it picks per-cell weights for each modality based on local prediction
  # strength, so cells where ATAC is informative get weighted toward ATAC.
  obj <- FindMultiModalNeighbors(
    obj,
    reduction.list = list("pca", names(obj@reductions)[
      grepl("lsi|svd", names(obj@reductions), ignore.case = TRUE)][1]),
    dims.list = list(1:30, 2:30),
    verbose = FALSE
  )
  graph_name <- "wsnn"
} else if (opt$signal == "rna") {
  obj <- prep_rna(obj)
  obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, verbose = FALSE)
  graph_name <- paste0(DefaultAssay(obj), "_snn")
} else if (opt$signal == "atac") {
  obj <- prep_atac(obj)
  atac_assay <- DefaultAssay(obj)
  obj <- FindNeighbors(obj,
                       reduction = grep("lsi|svd",
                                        names(obj@reductions),
                                        ignore.case = TRUE, value = TRUE)[1],
                       dims = 2:30, verbose = FALSE)
  graph_name <- paste0(atac_assay, "_snn")
} else {
  stop("unknown --signal: ", opt$signal)
}

# ------------------------------------------------------------------------
# Separate-assay safety check
# ------------------------------------------------------------------------
# If clustering on ATAC and no transferred-label column is present, the user
# probably skipped both the coembed-construct and Signac integrate_atac
# paths. Refuse with a pointer at coembed-construct (the simpler entry
# point) and the integrate_atac vignette. Skipped when the user has passed
# --use-existing-clusters (object is pre-coembedded).
if (!opt$`use-existing-clusters` && opt$signal == "atac") {
  if (!(opt$`transferred-label-col` %in% colnames(obj@meta.data))) {
    stop(sprintf(
      "--signal atac selected but column '%s' is missing from meta.data. ",
      opt$`transferred-label-col`),
      "If this is separate-assay scATAC that hasn't been integrated with ",
      "scRNA yet, run coembed-construct first to produce a coembed .rds, ",
      "then call this skill with --use-existing-clusters. Or, if you've ",
      "already done your own integrate_atac upstream, override with ",
      "--transferred-label-col <existing col>. ",
      "Vignette: https://stuartlab.org/signac/articles/integrate_atac")
  }
}

# ------------------------------------------------------------------------
# Resolution search + cluster assignment
# ------------------------------------------------------------------------
# Skipped when --use-existing-clusters: res_result was set above with the
# pre-clustered object and a sentinel trace.
if (!opt$`use-existing-clusters`) {
  res_result <- binary_search_resolution(
    obj, graph_name, n_cells, opt$`target-cluster-size`)
  obj <- res_result$obj
}

# Write the trace for audit.
writeLines(
  toJSON(list(signal = opt$signal,
              target_cluster_size = opt$`target-cluster-size`,
              chosen_resolution = res_result$resolution,
              n_cells = n_cells,
              iterations = res_result$trace),
         pretty = TRUE, auto_unbox = TRUE),
  file.path(output_dir, "resolution_trace.json")
)

# ------------------------------------------------------------------------
# Filter small clusters
# ------------------------------------------------------------------------
# For multiome / WNN, a cell is either in or out by cluster, so a single
# count check suffices. For separate-assay clustering on ATAC, we still get
# one cluster assignment per cell. The --min-cluster-cells spec says 20.

cluster_vec <- as.character(Idents(obj))
sizes <- table(cluster_vec)
keep_clusters <- names(sizes)[sizes >= opt$`min-cluster-cells`]
dropped_clusters <- setdiff(names(sizes), keep_clusters)
if (length(dropped_clusters) > 0) {
  message("[load_and_cluster] dropping ", length(dropped_clusters),
          " clusters below --min-cluster-cells=", opt$`min-cluster-cells`,
          ": ", paste(dropped_clusters, collapse = ","))
}

# ------------------------------------------------------------------------
# Metadata column selection
# ------------------------------------------------------------------------

auto_detect_metadata <- function(obj, cluster_col = "seurat_clusters") {
  md <- obj@meta.data
  keep <- c()
  for (col in colnames(md)) {
    if (col %in% c(cluster_col, "nCount_RNA", "nFeature_RNA",
                   "nCount_ATAC", "nFeature_ATAC", "percent.mt",
                   "orig.ident")) {
      # orig.ident sometimes useful but often noisy per-cell barcode source
      next
    }
    v <- md[[col]]
    if (is.numeric(v) && length(unique(v)) > 30) next  # continuous
    n_unique <- length(unique(v))
    if (n_unique >= 2 && n_unique <= 20) keep <- c(keep, col)
  }
  keep
}

if (!is.null(opt$`metadata-cols`)) {
  meta_cols <- strsplit(opt$`metadata-cols`, ",", fixed = TRUE)[[1]]
  meta_cols <- trimws(meta_cols)
  missing_cols <- setdiff(meta_cols, colnames(obj@meta.data))
  if (length(missing_cols) > 0) {
    stop("metadata columns not found in object: ",
         paste(missing_cols, collapse = ","))
  }
} else {
  meta_cols <- auto_detect_metadata(obj)
  message("[load_and_cluster] auto-detected metadata cols: ",
          paste(meta_cols, collapse = ","))
  if (!opt$yes && interactive()) {
    message("Pass --yes or --metadata-cols to accept non-interactively.")
  }
}

# ------------------------------------------------------------------------
# Write clusters.csv + per-cluster barcode files
# ------------------------------------------------------------------------

md_out <- obj@meta.data
md_out$barcode <- rownames(md_out)
md_out$seurat_cluster <- cluster_vec
# Always preserve `assay` if the object carries it — coembed outputs use
# this as the RNA/ATAC origin tag, and downstream aggregate_rna.R +
# call_peaks.R filter on it. Don't require the user to pass it in
# --metadata-cols just to keep the filter working.
preserve_cols <- unique(c("barcode", "seurat_cluster",
                          if ("assay" %in% colnames(md_out)) "assay",
                          meta_cols))
md_out <- md_out[md_out$seurat_cluster %in% keep_clusters,
                 preserve_cols, drop = FALSE]
write.csv(md_out, file.path(output_dir, "clusters.csv"), row.names = FALSE)

for (cl in keep_clusters) {
  cells <- md_out$barcode[md_out$seurat_cluster == cl]
  writeLines(cells,
             file.path(output_dir, "per_cluster_barcodes",
                       sprintf("cluster%s.txt", cl)))
}

# ------------------------------------------------------------------------
# Groups plan: cross-product of (cluster x meta_col_1 x meta_col_2 x ...)
# ------------------------------------------------------------------------
# When the user supplies multiple --metadata-cols (e.g. genotype,tissue),
# we want the FULL CROSS-PRODUCT within each cluster — not parallel
# per-column groups. The cross-product is the correct stratification when
# metadata axes interact (e.g. WT-spleen and WT-siiel have different
# chromatin landscapes that pooling would smooth out). The per-axis
# parallel pattern only makes sense for truly independent confounders,
# which is rarely the case in single-cell biology.
#
# Group spec format (groups_plan.json):
#   { name, cluster, metadata: {col1: val1, col2: val2, ...},
#     n_cells, cohort, group }
#
# `cohort` is the value of the col named in --cohort-col (default: first
# of --metadata-cols). It's the column-VALUE, not the column-NAME — that
# was a bug in the parallel-per-col implementation and would have produced
# wrong cohort labels in the build-taiji-input manifest.

groups <- list()

# Resolve cohort axis: either explicit --cohort-col, or the first
# --metadata-cols entry, or "all" if no metadata cols.
cohort_col <- opt$`cohort-col`
if (is.null(cohort_col) || cohort_col == "") {
  cohort_col <- if (length(meta_cols) > 0) meta_cols[1] else NA_character_
}
if (!is.na(cohort_col) && !(cohort_col %in% meta_cols)) {
  stop(sprintf(
    "--cohort-col '%s' is not among --metadata-cols (%s).",
    cohort_col, paste(meta_cols, collapse = ",")))
}

make_name <- function(cluster, metadata_pairs) {
  safe <- function(s) gsub("[^A-Za-z0-9._-]+", "_", as.character(s))
  if (length(metadata_pairs) == 0) {
    return(sprintf("cluster%s", safe(cluster)))
  }
  parts <- mapply(
    function(k, v) sprintf("%s-%s", safe(k), safe(v)),
    names(metadata_pairs), metadata_pairs, USE.NAMES = FALSE
  )
  sprintf("cluster%s__%s", safe(cluster), paste(parts, collapse = "__"))
}

# Build the (col -> sorted unique non-NA values) lookup once. We use this
# as input to expand.grid for the cross-product.
meta_value_lookup <- function(sub, cols) {
  out <- lapply(cols, function(col) {
    v <- unique(sub[[col]])
    sort(as.character(v[!is.na(v)]))
  })
  names(out) <- cols
  out
}

for (cl in keep_clusters) {
  sub <- md_out[md_out$seurat_cluster == cl, , drop = FALSE]

  if (length(meta_cols) == 0) {
    # No metadata stratification: one group per cluster.
    name <- make_name(cl, list())
    bcs <- sub$barcode
    writeLines(bcs,
               file.path(output_dir, "per_cluster_barcodes",
                         sprintf("%s.txt", name)))
    groups[[length(groups) + 1]] <- list(
      name = name, cluster = cl,
      metadata = NULL,
      n_cells = as.integer(nrow(sub)),
      cohort = "all", group = name
    )
    next
  }

  # Cross-product of all metadata cols' unique values within this cluster.
  vals <- meta_value_lookup(sub, meta_cols)
  if (any(sapply(vals, length) == 0)) {
    message("[load_and_cluster] cluster ", cl,
            ": one of the metadata cols has no non-NA values; skipping.")
    next
  }
  grid <- do.call(expand.grid,
                  c(vals, list(stringsAsFactors = FALSE)))

  for (i in seq_len(nrow(grid))) {
    combo <- as.list(grid[i, , drop = FALSE])
    # Ensure values are character (expand.grid sometimes coerces).
    combo <- lapply(combo, as.character)
    # Filter the per-cluster sub to rows matching every combo constraint.
    matching <- rep(TRUE, nrow(sub))
    for (col in names(combo)) {
      matching <- matching &
                  !is.na(sub[[col]]) &
                  as.character(sub[[col]]) == combo[[col]]
    }
    n <- sum(matching, na.rm = TRUE)
    if (n < opt$`min-cluster-cells`) next  # per-group floor

    name <- make_name(cl, combo)
    bcs <- sub$barcode[matching]
    writeLines(bcs,
               file.path(output_dir, "per_cluster_barcodes",
                         sprintf("%s.txt", name)))

    cohort_label <- if (!is.na(cohort_col)) {
      as.character(combo[[cohort_col]])
    } else {
      "all"
    }

    groups[[length(groups) + 1]] <- list(
      name = name, cluster = cl,
      metadata = combo,
      n_cells = as.integer(n),
      cohort = cohort_label,
      group = name
    )
  }
}

writeLines(
  toJSON(list(
    signal = opt$signal,
    chosen_resolution = res_result$resolution,
    kept_clusters = keep_clusters,
    dropped_clusters = dropped_clusters,
    metadata_cols = meta_cols,
    cohort_col = if (is.na(cohort_col)) NULL else cohort_col,
    min_cluster_cells = opt$`min-cluster-cells`,
    groups = groups
  ), pretty = TRUE, auto_unbox = TRUE),
  file.path(output_dir, "groups_plan.json")
)

message("[load_and_cluster] done. ",
        length(keep_clusters), " clusters kept, ",
        length(groups), " pseudobulk groups planned.")

# ------------------------------------------------------------------------
# QC UMAP — colored panels for clusters / assay / metadata
# ------------------------------------------------------------------------
# Generated unconditionally (default ON) once clustering has converged.
# The reduction is computed against the same neighbor graph used for
# clustering, so the UMAP visually reflects the partitions.
#
# Panels rendered:
#   - seurat_clusters (always)
#   - assay           (if a meta.data column 'assay' exists — the
#                      integrate_atac convention for co-embedded objects)
#   - each detected `--metadata-cols` value
#
# Cheap to compute (a few seconds for ≤100k cells); easy to skip via
# --no-plot if the user is iterating on a CI box without a graphics stack.
# Degrades gracefully when ggplot2 / patchwork aren't installed.

if (!opt$`no-plot`) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("[load_and_cluster] ggplot2 not available; skipping QC UMAP.")
  } else {
    qc_dir <- file.path(output_dir, "qc")
    dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

    # Compute UMAP using the same reduction we clustered on. When
    # --use-existing-clusters is set, prefer the object's existing UMAP
    # reduction so the QC plot matches what the user already has.
    obj_for_plot <- tryCatch({
      if (opt$`use-existing-clusters`) {
        existing_umap <- intersect(c("umap", "umap.wnn", "umap.atac", "umap.rna"),
                                   names(obj@reductions))[1]
        if (!is.na(existing_umap)) {
          obj   # already has a UMAP; no recompute
        } else if ("pca" %in% names(obj@reductions)) {
          RunUMAP(obj, reduction = "pca", dims = 1:30, verbose = FALSE)
        } else NULL
      } else if (opt$signal == "wnn") {
        RunUMAP(obj, nn.name = "weighted.nn",
                reduction.name = "umap.wnn", verbose = FALSE)
      } else if (opt$signal == "rna") {
        RunUMAP(obj, reduction = "pca", dims = 1:30, verbose = FALSE)
      } else if (opt$signal == "atac") {
        lsi_red <- grep("lsi|svd", names(obj@reductions),
                        ignore.case = TRUE, value = TRUE)[1]
        RunUMAP(obj, reduction = lsi_red, dims = 2:30,
                reduction.name = "umap.atac", verbose = FALSE)
      } else NULL
    }, error = function(e) {
      message("[load_and_cluster] UMAP failed (", conditionMessage(e),
              "); skipping QC plot.")
      NULL
    })

    if (!is.null(obj_for_plot)) {
      reduction_name <- if (opt$`use-existing-clusters`) {
        intersect(c("umap", "umap.wnn", "umap.atac", "umap.rna"),
                  names(obj_for_plot@reductions))[1]
      } else {
        switch(opt$signal,
               wnn  = "umap.wnn",
               atac = "umap.atac",
               "umap")
      }

      # Save the UMAP coords as CSV so downstream tools can re-render
      # without R (matplotlib, plotly, etc.). Cheap, ~few MB max.
      umap_emb <- Embeddings(obj_for_plot, reduction = reduction_name)
      umap_df <- data.frame(
        barcode = rownames(umap_emb),
        umap_1  = umap_emb[, 1],
        umap_2  = umap_emb[, 2],
        seurat_cluster = as.character(Idents(obj_for_plot)),
        stringsAsFactors = FALSE
      )
      write.csv(umap_df, file.path(qc_dir, "umap_coords.csv"), row.names = FALSE)

      # Pick group-by columns: clusters (always) + assay (if present)
      # + each user-chosen metadata col (only those that exist).
      groupby_cols <- c("seurat_clusters")
      if ("assay" %in% colnames(obj_for_plot@meta.data)) {
        groupby_cols <- c(groupby_cols, "assay")
      }
      groupby_cols <- c(groupby_cols,
                        intersect(meta_cols, colnames(obj_for_plot@meta.data)))
      groupby_cols <- unique(groupby_cols)

      panels <- list()
      for (col in groupby_cols) {
        p <- tryCatch(
          DimPlot(obj_for_plot, reduction = reduction_name, group.by = col,
                  label = (col == "seurat_clusters"), repel = TRUE) +
            ggplot2::ggtitle(col),
          error = function(e) {
            message("[load_and_cluster] DimPlot[", col, "] failed: ",
                    conditionMessage(e))
            NULL
          }
        )
        if (!is.null(p)) panels[[length(panels) + 1]] <- p
      }

      if (length(panels) > 0) {
        out_png <- file.path(qc_dir, "umap.png")
        if (requireNamespace("patchwork", quietly = TRUE) &&
            length(panels) > 1) {
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
        message("[load_and_cluster] QC UMAP -> ", out_png)
      }
    }
  }
}
