#!/usr/bin/env Rscript
# build_atac_from_fragments.R
# Build and QC-filter a ChromatinAssay Seurat object from a raw fragments.tsv.gz
# when no pre-built .rds exists.
#
# Pipeline:
#   1. CountFragments()       per-barcode fragment tally; apply --min-fragments threshold
#   2. Pre-filter fragments   zcat|awk to a cell-specific BED for retained barcodes
#   3. MACS3                  call peaks on the filtered BED
#   4. FeatureMatrix()        quantify peaks x cells
#   5. CreateChromatinAssay() + CreateSeuratObject()
#   6. NucleosomeSignal()     always (needs only fragments)
#   7. TSSEnrichment()        if --genome supplied (needs EnsDb)
#   8. blacklist_ratio        if --genome supplied (needs Signac blacklist data)
#   9. Apply QC filters; emit filtered .rds, qc_summary.json, violin plots

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(ggplot2)
  library(jsonlite)
})

option_list <- list(
  make_option("--fragments",      type = "character",
              help = "Path to fragments.tsv.gz (must have .tbi sibling index)."),
  make_option("--output",         type = "character",
              help = "Output .rds path for the QC-filtered Seurat object."),
  make_option("--genome",         type = "character", default = NULL,
              help = paste0("Genome build (hg19, hg38, mm10, mm39). ",
                            "Enables TSS enrichment and blacklist-ratio QC. ",
                            "Requires the matching EnsDb.* package and Signac blacklist data.")),
  make_option("--min-fragments",  type = "integer",   default = 1000L,
              help = paste0("Min fragments per barcode (passed_filters from CountFragments). ",
                            "Barcodes below this threshold are treated as empty. Default: 1000.")),
  make_option("--max-fragments",  type = "integer",   default = NULL,
              help = "Max fragments per barcode; NULL = no upper limit (default)."),
  make_option("--min-count",      type = "integer",   default = 1000L,
              help = "nCount_ATAC min after FeatureMatrix quantification. Default: 1000."),
  make_option("--min-features",   type = "integer",   default = 500L,
              help = "nFeature_ATAC (# peaks) min after FeatureMatrix. Default: 500."),
  make_option("--tss-min",        type = "double",    default = 2.0,
              help = "TSS.enrichment minimum (only if --genome supplied). Default: 2.0."),
  make_option("--nucleosome-max", type = "double",    default = 4.0,
              help = "nucleosome_signal maximum. Default: 4.0."),
  make_option("--blacklist-max",  type = "double",    default = 0.05,
              help = "blacklist_ratio maximum (only if --genome supplied). Default: 0.05."),
  make_option("--macs-path",      type = "character", default = NULL,
              help = "Path to macs2 or macs3 binary. Auto-detected from PATH if omitted."),
  make_option("--assay-name",     type = "character", default = "ATAC",
              help = "ChromatinAssay name in the output Seurat object. Default: ATAC."),
  make_option("--group-by",       type = "character", default = NULL,
              help = "meta.data column for violin-plot grouping (optional).")
)
opt <- parse_args(OptionParser(option_list = option_list))

for (k in c("fragments", "output")) {
  if (is.null(opt[[k]])) stop("missing required --", k)
}
if (!file.exists(opt$fragments)) stop("fragments file not found: ", opt$fragments)
tbi_path <- paste0(opt$fragments, ".tbi")
if (!file.exists(tbi_path)) {
  stop(".tbi tabix index not found: ", tbi_path,
       "\nRun: tabix -p bed ", opt$fragments)
}

# Auto-append .rds if extension is missing.
if (!grepl("\\.rds$", opt$output, ignore.case = TRUE)) {
  opt$output <- paste0(opt$output, ".rds")
  message("[build_atac] NOTE: .rds appended to --output path -> ", opt$output)
}

output_dir <- dirname(opt$output)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
qc_dir <- file.path(output_dir, "qc")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# MACS auto-detection
# ---------------------------------------------------------------------------

resolve_macs <- function(user_path) {
  if (!is.null(user_path) && nchar(user_path) > 0) {
    if (!file.exists(user_path)) stop("--macs-path not found: ", user_path)
    return(user_path)
  }
  for (bin in c("macs3", "macs2")) {
    found <- Sys.which(bin)
    if (nchar(found) > 0) {
      message("[build_atac] macs binary: ", found)
      return(found)
    }
  }
  stop("macs2/macs3 not found on PATH. Install MACS3 or supply --macs-path.")
}

macs_path <- resolve_macs(opt$`macs-path`)

# ---------------------------------------------------------------------------
# Effective genome size (ENCODE / Signac defaults)
# ---------------------------------------------------------------------------

effective_genome_size <- function(genome) {
  sizes <- list(hg19 = 2.7e9, hg38 = 2.9e9, mm9 = 1.87e9, mm10 = 1.87e9, mm39 = 1.87e9)
  sz <- sizes[[tolower(genome)]]
  if (is.null(sz)) {
    stop("unsupported --genome '", genome,
         "'. Known: ", paste(names(sizes), collapse = ", "),
         ". Add the effective genome size for your build.")
  }
  sz
}

# ---------------------------------------------------------------------------
# Step 1 — CountFragments: per-barcode tally
# ---------------------------------------------------------------------------

message("[build_atac] Step 1: CountFragments() on ", basename(opt$fragments))
frag_counts <- CountFragments(fragments = opt$fragments)
message(sprintf("[build_atac]   total barcodes in fragment file: %d",
                nrow(frag_counts)))

# passed_filters is the deduplicated-read count; frequency_count includes
# duplicates. Use passed_filters as the quality metric, falling back to
# frequency_count if the column is absent (older Signac versions).
count_col <- if ("passed_filters" %in% colnames(frag_counts)) "passed_filters" else "frequency_count"
message("[build_atac]   using column '", count_col, "' for barcode filtering")

keep_bc <- frag_counts$CB[frag_counts[[count_col]] >= opt$`min-fragments`]
if (!is.null(opt$`max-fragments`)) {
  keep_bc <- keep_bc[frag_counts[[count_col]][frag_counts$CB %in% keep_bc] <= opt$`max-fragments`]
}
n_total  <- nrow(frag_counts)
n_kept   <- length(keep_bc)
n_empty  <- n_total - n_kept
message(sprintf("[build_atac]   barcodes >= %d fragments: %d / %d  (%.1f%% retained)",
                opt$`min-fragments`, n_kept, n_total, 100 * n_kept / n_total))
message(sprintf("[build_atac]   discarded as empty: %d (%.1f%%)",
                n_empty, 100 * n_empty / n_total))
if (n_kept == 0) stop("No barcodes passed --min-fragments=", opt$`min-fragments`,
                      ". Lower the threshold or check the fragments file.")

# Save barcode rank plot
rank_df <- frag_counts[order(frag_counts[[count_col]], decreasing = TRUE), ]
rank_df$rank <- seq_len(nrow(rank_df))
p_rank <- ggplot2::ggplot(rank_df, ggplot2::aes(x = rank, y = .data[[count_col]])) +
  ggplot2::geom_line(colour = "#1f77b4") +
  ggplot2::geom_hline(yintercept = opt$`min-fragments`, linetype = "dashed",
                      colour = "#d62728", linewidth = 0.5) +
  ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
  ggplot2::labs(title = "Barcode-rank plot (log-log)",
                x = "Barcode rank", y = count_col) +
  ggplot2::theme_bw(base_size = 11)
ggplot2::ggsave(file.path(qc_dir, "barcode_rank_plot.png"),
                plot = p_rank, width = 7, height = 5, dpi = 120, bg = "white")
message("[build_atac]   -> qc/barcode_rank_plot.png")

# ---------------------------------------------------------------------------
# Step 2 — Pre-filter fragments to retained barcodes (zcat|awk)
# ---------------------------------------------------------------------------

message("[build_atac] Step 2: extracting fragments for ", n_kept, " retained barcodes ...")
bc_file <- tempfile(fileext = ".txt")
writeLines(keep_bc, bc_file)
tmp_bed <- tempfile(pattern = "atac_frags_", fileext = ".bed")
on.exit(unlink(c(bc_file, tmp_bed), force = TRUE), add = TRUE)

frag_cmd <- sprintf(
  "zcat %s | awk -v bcfile=%s 'BEGIN{while((getline l<bcfile)>0) keep[l]=1} /^#/{next} $4 in keep {print $1\"\\t\"$2\"\\t\"$3}' > %s",
  shQuote(opt$fragments), shQuote(bc_file), shQuote(tmp_bed)
)
rc <- system(frag_cmd)
if (rc != 0 || !file.exists(tmp_bed)) {
  stop("fragment pre-filtering failed (exit code ", rc,
       "). Is the file bgzipped and is zcat/awk on PATH?")
}
n_frags_extracted <- tryCatch(
  as.integer(system(sprintf("wc -l < %s", shQuote(tmp_bed)), intern = TRUE)),
  error = function(e) NA_integer_
)
message(sprintf("[build_atac]   extracted %s fragment rows",
                format(n_frags_extracted, big.mark = ",")))

# ---------------------------------------------------------------------------
# Step 3 — CallPeaks: MACS3 on the retained-barcode BED
# ---------------------------------------------------------------------------

message("[build_atac] Step 3: calling peaks via MACS3 ...")
gsize <- if (!is.null(opt$genome)) effective_genome_size(opt$genome) else 2.7e9
if (is.null(opt$genome)) message("[build_atac]   WARN: --genome not supplied; using gsize=2.7e9 (hg19). ",
                                  "Supply --genome for the correct effective genome size.")

tmp_macs <- tempfile(pattern = "macs_out_")
dir.create(tmp_macs, showWarnings = FALSE)
on.exit(unlink(tmp_macs, recursive = TRUE), add = TRUE)

macs_cmd <- sprintf(
  "%s callpeak -t %s -g %.2e -f BED --nomodel --extsize 200 --shift -100 -n peaks --outdir %s 2>&1",
  shQuote(macs_path), shQuote(tmp_bed), gsize, shQuote(tmp_macs)
)
rc_macs <- system(macs_cmd)
np_file <- file.path(tmp_macs, "peaks_peaks.narrowPeak")
if (rc_macs != 0 || !file.exists(np_file) || file.size(np_file) == 0) {
  stop("MACS3 peak calling failed (exit code ", rc_macs,
       "). Check MACS3 installation and fragment BED content.")
}

np_df <- read.table(np_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(np_df) <- c("chrom", "chromStart", "chromEnd", "name", "score",
                     "strand", "signalValue", "pValue", "qValue", "peak")
n_peaks <- nrow(np_df)
message(sprintf("[build_atac]   %d peaks called", n_peaks))

peaks_gr <- GenomicRanges::GRanges(
  seqnames = np_df$chrom,
  ranges   = IRanges::IRanges(np_df$chromStart + 1L, np_df$chromEnd),
  name     = np_df$name,
  score    = np_df$score,
  fold_change              = np_df$signalValue,
  neg_log10pvalue_summit   = np_df$pValue,
  neg_log10qvalue_summit   = np_df$qValue,
  relative_summit_position = np_df$peak
)

# ---------------------------------------------------------------------------
# Step 4 — FeatureMatrix: quantify peaks x cells
# ---------------------------------------------------------------------------

message("[build_atac] Step 4: FeatureMatrix() quantifying ", n_peaks,
        " peaks x ", n_kept, " cells ...")
frags_obj <- tryCatch(
  CreateFragmentObject(path = opt$fragments, cells = keep_bc, validate.fragments = FALSE),
  error = function(e) stop("CreateFragmentObject failed: ", conditionMessage(e),
                            "\nCheck that the .tbi index is current.")
)

counts_mat <- FeatureMatrix(
  fragments = frags_obj,
  features  = peaks_gr,
  cells     = keep_bc
)
message(sprintf("[build_atac]   count matrix: %d peaks x %d cells",
                nrow(counts_mat), ncol(counts_mat)))

# ---------------------------------------------------------------------------
# Step 5 — CreateChromatinAssay + CreateSeuratObject
# ---------------------------------------------------------------------------

message("[build_atac] Step 5: CreateChromatinAssay() + CreateSeuratObject()")
chrom_assay <- CreateChromatinAssay(
  counts    = counts_mat,
  sep       = c(":", "-"),
  fragments = opt$fragments,
  min.cells = 1,
  min.features = 1
)

obj <- CreateSeuratObject(counts = chrom_assay, assay = opt$`assay-name`)
message(sprintf("[build_atac]   Seurat object: %d cells", ncol(obj)))
n_before <- ncol(obj)

# ---------------------------------------------------------------------------
# Step 6 — NucleosomeSignal (always; needs only fragments)
# ---------------------------------------------------------------------------

message("[build_atac] Step 6: NucleosomeSignal() ...")
obj <- tryCatch(
  NucleosomeSignal(obj),
  error = function(e) {
    message("[build_atac]   WARN: NucleosomeSignal failed: ", conditionMessage(e),
            "\n  nucleosome_signal filter will be skipped.")
    obj
  }
)

# ---------------------------------------------------------------------------
# Step 7 — TSSEnrichment and blacklist_ratio (requires --genome)
# ---------------------------------------------------------------------------

if (!is.null(opt$genome)) {
  genome_lc <- tolower(opt$genome)

  # EnsDb dispatch
  ensdb_pkg <- list(
    hg38 = "EnsDb.Hsapiens.v86",
    hg19 = "EnsDb.Hsapiens.v86",
    mm10 = "EnsDb.Mmusculus.v79",
    mm39 = "EnsDb.Mmusculus.v79"
  )[[genome_lc]]
  if (is.null(ensdb_pkg)) {
    message("[build_atac]   WARN: no EnsDb mapping for --genome '", opt$genome,
            "'. TSS enrichment and blacklist ratio skipped.")
  } else if (!requireNamespace(ensdb_pkg, quietly = TRUE)) {
    message("[build_atac]   WARN: package '", ensdb_pkg, "' not installed. ",
            "TSS enrichment and blacklist ratio skipped.")
  } else {
    message("[build_atac] Step 7: TSSEnrichment() and blacklist_ratio (genome=", opt$genome, ") ...")
    ensdb <- get(ensdb_pkg, envir = asNamespace(ensdb_pkg))
    annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
    # seqlevelsStyle<- may attempt a UCSC network fetch on HPC nodes without internet.
    tryCatch(
      seqlevelsStyle(annotations) <- "UCSC",
      error = function(e) {
        message("[build_atac]   WARN: seqlevelsStyle UCSC fetch failed (",
                conditionMessage(e), "). Applying manual chr-prefix fallback.")
        lvls <- seqlevels(annotations)
        new_lvls <- ifelse(grepl("^chr", lvls), lvls,
                           ifelse(lvls == "MT", "chrM", paste0("chr", lvls)))
        seqlevels(annotations) <<- new_lvls
      }
    )
    genome(annotations) <- opt$genome
    Annotation(obj) <- annotations

    obj <- tryCatch(
      TSSEnrichment(obj, fast = FALSE),
      error = function(e) {
        message("[build_atac]   WARN: TSSEnrichment failed: ", conditionMessage(e))
        obj
      }
    )

    # Blacklist ratio via Signac's built-in blacklist data
    bl_name <- paste0("blacklist_", genome_lc)
    bl_obj  <- tryCatch(
      get(bl_name, envir = asNamespace("Signac")),
      error = function(e) NULL
    )
    if (!is.null(bl_obj)) {
      obj$blacklist_ratio <- FractionCountsInRegion(
        object  = obj,
        assay   = opt$`assay-name`,
        regions = bl_obj
      )
    } else {
      message("[build_atac]   WARN: Signac::", bl_name, " not found. ",
              "blacklist_ratio filter skipped.")
    }
  }
} else {
  message("[build_atac] Step 7: --genome not supplied; skipping TSS enrichment and blacklist ratio.")
}

# ---------------------------------------------------------------------------
# Violin plots — BEFORE filtering
# ---------------------------------------------------------------------------

atac_qc_cols <- c(
  paste0("nCount_",   opt$`assay-name`),
  paste0("nFeature_", opt$`assay-name`),
  "nucleosome_signal", "TSS.enrichment", "blacklist_ratio"
)
present_cols <- intersect(atac_qc_cols, colnames(obj@meta.data))

grp <- opt$`group-by`
if (!is.null(grp) && !(grp %in% colnames(obj@meta.data))) grp <- NULL

plot_atac_violins <- function(obj, title_suffix, grp, thresholds) {
  present <- intersect(atac_qc_cols, colnames(obj@meta.data))
  lapply(present, function(col) {
    p <- VlnPlot(obj, features = col, group.by = grp, pt.size = 0) +
      ggplot2::ggtitle(paste0(col, " ", title_suffix)) +
      ggplot2::theme(legend.position = "none",
                     axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    if (!is.null(thresholds[[col]])) {
      p <- p + ggplot2::geom_hline(
        yintercept = thresholds[[col]]$value, linetype = "dashed",
        colour = "#d62728", linewidth = 0.5)
    }
    p
  })
}

ncount_col   <- paste0("nCount_",   opt$`assay-name`)
nfeature_col <- paste0("nFeature_", opt$`assay-name`)
thresholds <- list(
  nucleosome_signal = list(value = opt$`nucleosome-max`),
  TSS.enrichment    = list(value = opt$`tss-min`),
  blacklist_ratio   = list(value = opt$`blacklist-max`)
)
thresholds[[ncount_col]]   <- list(value = opt$`min-count`)
thresholds[[nfeature_col]] <- list(value = opt$`min-features`)

message("[build_atac] plotting QC violins before filtering ...")
before_plots <- plot_atac_violins(obj, "(before)", grp, thresholds)
if (length(before_plots) > 0 && requireNamespace("patchwork", quietly = TRUE)) {
  ggplot2::ggsave(
    file.path(qc_dir, "violin_before.png"),
    plot  = patchwork::wrap_plots(before_plots, ncol = min(3L, length(before_plots))),
    width = max(7, 5 * min(3L, length(before_plots))), height = 5,
    dpi   = 120, bg = "white"
  )
  message("[build_atac]   -> qc/violin_before.png")
}

# ---------------------------------------------------------------------------
# Step 9 — Apply QC filters
# ---------------------------------------------------------------------------

message("[build_atac] Step 9: applying QC filters ...")

md <- obj@meta.data
keep <- rep(TRUE, ncol(obj))

n_fail <- list()

if (ncount_col %in% colnames(md)) {
  fail <- md[[ncount_col]] < opt$`min-count`
  n_fail$ncount_lo <- sum(fail)
  keep <- keep & !fail
  message(sprintf("  %s < %d : %d cells", ncount_col, opt$`min-count`, n_fail$ncount_lo))
}

if (nfeature_col %in% colnames(md)) {
  fail <- md[[nfeature_col]] < opt$`min-features`
  n_fail$nfeature_lo <- sum(fail)
  keep <- keep & !fail
  message(sprintf("  %s < %d : %d cells", nfeature_col, opt$`min-features`, n_fail$nfeature_lo))
}

if ("nucleosome_signal" %in% colnames(md)) {
  fail <- md$nucleosome_signal >= opt$`nucleosome-max`
  n_fail$nucleosome_hi <- sum(fail)
  keep <- keep & !fail
  message(sprintf("  nucleosome_signal >= %.1f : %d cells", opt$`nucleosome-max`, n_fail$nucleosome_hi))
} else {
  message("  nucleosome_signal: not in meta.data — skipped")
}

if ("TSS.enrichment" %in% colnames(md)) {
  fail <- md$TSS.enrichment < opt$`tss-min`
  n_fail$tss_lo <- sum(fail)
  keep <- keep & !fail
  message(sprintf("  TSS.enrichment < %.1f : %d cells", opt$`tss-min`, n_fail$tss_lo))
} else {
  message("  TSS.enrichment: not in meta.data — skipped (supply --genome to enable)")
}

if ("blacklist_ratio" %in% colnames(md)) {
  fail <- md$blacklist_ratio >= opt$`blacklist-max`
  n_fail$blacklist_hi <- sum(fail)
  keep <- keep & !fail
  message(sprintf("  blacklist_ratio >= %.3f : %d cells", opt$`blacklist-max`, n_fail$blacklist_hi))
} else {
  message("  blacklist_ratio: not in meta.data — skipped (supply --genome to enable)")
}

n_after   <- sum(keep)
n_removed <- n_before - n_after
message(sprintf("  total removed (any filter): %d / %d  (%.1f%%)",
                n_removed, n_before, 100 * n_removed / n_before))
message(sprintf("  retained: %d cells", n_after))

if (n_after == 0) stop("All cells removed by QC filters. Loosen thresholds.")

obj <- subset(obj, cells = colnames(obj)[keep])

# ---------------------------------------------------------------------------
# Violin plots — AFTER filtering
# ---------------------------------------------------------------------------

message("[build_atac] plotting QC violins after filtering ...")
after_plots <- plot_atac_violins(obj, "(after)", grp, thresholds)
if (length(after_plots) > 0 && requireNamespace("patchwork", quietly = TRUE)) {
  ggplot2::ggsave(
    file.path(qc_dir, "violin_after.png"),
    plot  = patchwork::wrap_plots(after_plots, ncol = min(3L, length(after_plots))),
    width = max(7, 5 * min(3L, length(after_plots))), height = 5,
    dpi   = 120, bg = "white"
  )
  message("[build_atac]   -> qc/violin_after.png")

  all_panels <- c(before_plots, after_plots)
  ggplot2::ggsave(
    file.path(qc_dir, "violin_comparison.png"),
    plot  = patchwork::wrap_plots(all_panels, ncol = min(3L, length(before_plots))),
    width = max(7, 5 * min(3L, length(before_plots))), height = 10,
    dpi   = 120, bg = "white"
  )
  message("[build_atac]   -> qc/violin_comparison.png")
}

# ---------------------------------------------------------------------------
# Save + summary
# ---------------------------------------------------------------------------

saveRDS(obj, opt$output)
message("[build_atac] saved -> ", opt$output)

summary_out <- list(
  fragments_input   = normalizePath(opt$fragments, mustWork = FALSE),
  output            = normalizePath(opt$output,    mustWork = FALSE),
  genome            = opt$genome,
  n_barcodes_total  = n_total,
  n_barcodes_kept_step1 = n_kept,
  n_peaks_called    = n_peaks,
  n_cells_before_qc = n_before,
  n_cells_after_qc  = n_after,
  n_removed         = n_removed,
  pct_removed       = round(100 * n_removed / n_before, 2),
  filters = list(
    min_fragments    = opt$`min-fragments`,
    max_fragments    = opt$`max-fragments`,
    min_count        = opt$`min-count`,
    min_features     = opt$`min-features`,
    nucleosome_max   = opt$`nucleosome-max`,
    tss_min          = opt$`tss-min`,
    blacklist_max    = opt$`blacklist-max`
  ),
  n_fail = n_fail,
  metrics_computed = list(
    nucleosome_signal = "nucleosome_signal" %in% colnames(obj@meta.data),
    tss_enrichment    = "TSS.enrichment"    %in% colnames(obj@meta.data),
    blacklist_ratio   = "blacklist_ratio"   %in% colnames(obj@meta.data)
  )
)
writeLines(
  toJSON(summary_out, pretty = TRUE, auto_unbox = TRUE),
  file.path(output_dir, "atac_qc_summary.json")
)
message("[build_atac] summary -> atac_qc_summary.json")
message("[build_atac] done.")
