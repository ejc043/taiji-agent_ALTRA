#!/usr/bin/env Rscript
# filter_rna.R — cell-level QC filtering for scRNA-seq Seurat objects.
#
# Filters: nFeature_RNA (min,max) strictly, percent.mt < max.
# Matches the lab reference (coembed_preprocess.R): nFeature > 200 & < 5000 & percent.mt < 10.
# nCount_RNA filter is intentionally absent from the reference and disabled by default here.
# Outputs: filtered .rds, qc_summary.json, violin plots before and after.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(jsonlite)
})

option_list <- list(
  make_option("--input",             type = "character", help = "Input scRNA Seurat .rds"),
  make_option("--output",            type = "character", help = "Filtered .rds output path"),
  make_option("--min-features",      type = "integer",   default = 200L,
              help = "nFeature_RNA strictly > this (exclusive, matching lab reference)"),
  make_option("--max-features",      type = "integer",   default = 5000L,
              help = "nFeature_RNA strictly < this (exclusive, matching lab reference)"),
  make_option("--max-pct-mt",        type = "double",    default = 10.0),
  make_option("--mt-pattern",        type = "character", default = "^mt-",
              help = "Regex for mitochondrial genes (default ^mt- for mouse; use ^MT- for human)"),
  make_option("--group-by",          type = "character", default = NULL,
              help = "Metadata column for violin plot grouping (e.g. sample_id, tissue)")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (k in c("input", "output")) {
  if (is.null(opt[[k]])) stop("missing required --", k)
}

output_dir <- dirname(opt$output)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
qc_dir <- file.path(output_dir, "qc")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

message("[sc-qc] loading ", opt$input)
obj <- readRDS(opt$input)
n_before <- ncol(obj)
message(sprintf("[sc-qc]   %d cells before filtering", n_before))

# Compute percent.mt if absent
if (!"percent.mt" %in% colnames(obj@meta.data)) {
  message("[sc-qc]   computing percent.mt with pattern '", opt$`mt-pattern`, "'")
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = opt$`mt-pattern`)
}

# -------------------------------------------------------------------------
# Violin plots — BEFORE
# -------------------------------------------------------------------------
plot_violins <- function(obj, title_suffix, group_by = NULL) {
  cols <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  cols <- intersect(cols, colnames(obj@meta.data))
  plots <- lapply(cols, function(col) {
    p <- VlnPlot(obj, features = col, group.by = group_by,
                 pt.size = 0, log = (col == "nCount_RNA")) +
      ggplot2::ggtitle(paste0(col, " ", title_suffix)) +
      ggplot2::theme(legend.position = "none",
                     axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    if (col == "nFeature_RNA") {
      p <- p +
        ggplot2::geom_hline(yintercept = opt$`min-features`, linetype = "dashed",
                            colour = "#d62728", linewidth = 0.5) +
        ggplot2::geom_hline(yintercept = opt$`max-features`, linetype = "dashed",
                            colour = "#d62728", linewidth = 0.5)
    } else if (col == "percent.mt") {
      p <- p +
        ggplot2::geom_hline(yintercept = opt$`max-pct-mt`, linetype = "dashed",
                            colour = "#d62728", linewidth = 0.5)
    }
    p
  })
  plots
}

grp <- if (!is.null(opt$`group-by`) && opt$`group-by` %in% colnames(obj@meta.data)) {
  opt$`group-by`
} else if ("sample_id" %in% colnames(obj@meta.data)) {
  "sample_id"
} else if ("tissue" %in% colnames(obj@meta.data)) {
  "tissue"
} else {
  NULL
}

message("[sc-qc] plotting QC violins before filtering (group_by=", grp %||% "NULL", ")")
before_plots <- plot_violins(obj, "(before)", group_by = grp)
if (requireNamespace("patchwork", quietly = TRUE)) {
  ggplot2::ggsave(
    file.path(qc_dir, "violin_before.png"),
    plot = patchwork::wrap_plots(before_plots, ncol = 3),
    width = 14, height = 5, dpi = 120, bg = "white"
  )
  message("[sc-qc]   -> qc/violin_before.png")
}

# -------------------------------------------------------------------------
# Apply filters
# -------------------------------------------------------------------------
# Strict inequalities match the lab reference: nFeature > min & nFeature < max & percent.mt < max
# nCount_RNA is intentionally not filtered (absent from the reference pipeline).
n_fail_feat_lo <- sum(obj$nFeature_RNA <= opt$`min-features`)
n_fail_feat_hi <- sum(obj$nFeature_RNA >= opt$`max-features`)
n_fail_mt      <- sum(obj$percent.mt   >= opt$`max-pct-mt`)

keep <-
  obj$nFeature_RNA > opt$`min-features` &
  obj$nFeature_RNA < opt$`max-features` &
  obj$percent.mt   < opt$`max-pct-mt`

n_after <- sum(keep)
n_removed <- n_before - n_after

n_after   <- sum(keep)
n_removed <- n_before - n_after

message(sprintf("[sc-qc] filter summary (nFeature > %d & < %d & percent.mt < %.0f%%):",
                opt$`min-features`, opt$`max-features`, opt$`max-pct-mt`))
message(sprintf("  nFeature_RNA <= %d : %d cells", opt$`min-features`, n_fail_feat_lo))
message(sprintf("  nFeature_RNA >= %d : %d cells", opt$`max-features`, n_fail_feat_hi))
message(sprintf("  percent.mt   >= %.0f%%: %d cells", opt$`max-pct-mt`, n_fail_mt))
message(sprintf("  total removed (any filter): %d / %d  (%.1f%%)",
                n_removed, n_before, 100 * n_removed / n_before))
message(sprintf("  retained: %d cells", n_after))

obj <- subset(obj, cells = colnames(obj)[keep])

# -------------------------------------------------------------------------
# Violin plots — AFTER
# -------------------------------------------------------------------------
message("[sc-qc] plotting QC violins after filtering")
after_plots <- plot_violins(obj, "(after)", group_by = grp)
if (requireNamespace("patchwork", quietly = TRUE)) {
  ggplot2::ggsave(
    file.path(qc_dir, "violin_after.png"),
    plot = patchwork::wrap_plots(after_plots, ncol = 3),
    width = 14, height = 5, dpi = 120, bg = "white"
  )
  message("[sc-qc]   -> qc/violin_after.png")
}

# side-by-side before/after for easy comparison
all_panels <- c(before_plots, after_plots)
if (requireNamespace("patchwork", quietly = TRUE) && length(all_panels) > 0) {
  ggplot2::ggsave(
    file.path(qc_dir, "violin_comparison.png"),
    plot = patchwork::wrap_plots(all_panels, ncol = 3),
    width = 14, height = 10, dpi = 120, bg = "white"
  )
  message("[sc-qc]   -> qc/violin_comparison.png")
}

# -------------------------------------------------------------------------
# Save + summary
# -------------------------------------------------------------------------
saveRDS(obj, opt$output)
message("[sc-qc] saved -> ", opt$output)

summary_out <- list(
  input          = normalizePath(opt$input, mustWork = FALSE),
  output         = normalizePath(opt$output, mustWork = FALSE),
  n_cells_before = n_before,
  n_cells_after  = n_after,
  n_removed      = n_removed,
  pct_removed    = round(100 * n_removed / n_before, 2),
  filters = list(
    min_features = opt$`min-features`,
    max_features = opt$`max-features`,
    max_pct_mt   = opt$`max-pct-mt`,
    note         = "nCount_RNA not filtered (matches lab reference coembed_preprocess.R)"
  ),
  n_fail = list(
    feat_lo = n_fail_feat_lo,
    feat_hi = n_fail_feat_hi,
    mt      = n_fail_mt
  )
)
writeLines(
  toJSON(summary_out, pretty = TRUE, auto_unbox = TRUE),
  file.path(output_dir, "rna_qc_summary.json")
)
message("[sc-qc] summary -> rna_qc_summary.json")
message("[sc-qc] done.")
