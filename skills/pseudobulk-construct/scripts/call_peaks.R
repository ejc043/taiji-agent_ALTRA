#!/usr/bin/env Rscript
# call_peaks.R
#
# Per-group ATAC peak calling via Signac::CallPeaks (which itself shells out
# to macs2/macs3). Replaces the prior gunzip|awk|macs2 Python pipe so we get
# Signac's canonical fragment-object plumbing for free, and so peak calls
# survive the common barcode-suffix mismatch between Seurat objects and the
# raw 10x fragments file.
#
# Per group (one row in groups_plan.json):
#   1. Filter the cell list to the cluster (+ optional metadata col/value).
#   2. If a meta.data column 'assay' exists, restrict to assay == 'ATAC' —
#      this matches the integrate_atac convention where one co-embedded
#      object holds both modalities.
#   3. Reconcile barcodes against the fragments file. Seurat objects often
#      carry suffixes like 'AAACGCT-1' or 'ATAC_AAACGCT' that 10x's
#      fragments.tsv.gz does NOT have. We sample unique barcodes from the
#      fragments file once at startup, then for each group try the candidate
#      strip patterns and pick the one with highest overlap. Fail loudly
#      (≥95% overlap required) so the user knows when fragments and the
#      object can't be matched.
#   4. RenameCells in the SUBSETTED object (avoids creating duplicate names
#      when both RNA and ATAC sides of a co-embedded object share bare
#      barcodes), attach a fragment-restricted CreateFragmentObject, and call
#      Signac::CallPeaks.
#   5. Write <group>_peaks.rds and <group>.narrowPeak (ENCODE 10-col format).
#
# Idempotent: a group whose <group>_peaks.rds already exists is skipped, so
# partial reruns don't redo expensive macs invocations.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(jsonlite)
})

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

option_list <- list(
  make_option("--input",        type = "character",
              help = "Single-cell object (.rds/.h5ad); same one used by load_and_cluster.R."),
  make_option("--clusters",     type = "character",
              help = "clusters.csv produced by load_and_cluster.R."),
  make_option("--groups",       type = "character",
              help = "groups_plan.json produced by load_and_cluster.R."),
  make_option("--fragments",    type = "character",
              help = "Path to fragments.tsv.gz (must be tabix-indexed; .tbi sibling required by Signac)."),
  make_option("--output-dir",   type = "character",
              help = "Where to write <group>_peaks.rds and <group>.narrowPeak."),
  make_option("--genome",       type = "character",
              help = "Genome tag (hg19, hg38, mm9, mm10, mm39)."),
  make_option("--macs-path",    type = "character",
              help = "Full path to macs2 or macs3 binary."),
  make_option("--assay-col",    type = "character", default = "assay",
              help = "meta.data column flagging modality in co-embedded objects (value 'ATAC' is the ATAC subset). Default: 'assay'."),
  make_option("--min-overlap",  type = "double", default = 0.95,
              help = "Minimum fraction of group cells that must match the fragments file after suffix-stripping. Default: 0.95."),
  make_option("--frag-sample",  type = "integer", default = 5000000L,
              help = "Number of fragment rows to scan when collecting unique barcodes for the overlap check. Default: 5,000,000.")
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("input", "clusters", "groups", "fragments", "output-dir",
              "genome", "macs-path")
for (k in required) {
  if (is.null(opt[[k]])) stop("missing required --", k)
}
if (!file.exists(opt$fragments)) {
  stop("fragments file not found: ", opt$fragments)
}
dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Object loader (kept in lock-step with the loaders in the sibling R scripts;
# duplicated to avoid Rscript source() path gymnastics).
# ---------------------------------------------------------------------------

load_input <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    obj <- readRDS(path)
    if (!inherits(obj, "Seurat")) stop("not a Seurat object")
    return(obj)
  }
  if (ext == "h5ad") {
    if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
      stop("SeuratDisk required for .h5ad input — install manually if needed.")
    }
    h5seurat <- sub("\\.h5ad$", ".h5seurat", path, ignore.case = TRUE)
    if (!file.exists(h5seurat)) {
      SeuratDisk::Convert(path, dest = "h5seurat", overwrite = FALSE)
    }
    return(SeuratDisk::LoadH5Seurat(h5seurat))
  }
  stop("unsupported extension: .", ext)
}

# ---------------------------------------------------------------------------
# Genome size (Signac::CallPeaks needs an effective.genome.size; numbers from
# the ENCODE pipeline / Signac docs).
# ---------------------------------------------------------------------------

effective_genome_size <- function(genome) {
  table <- list(
    hg19 = 2.7e9,
    hg38 = 2.9e9,
    mm9  = 1.87e9,
    mm10 = 1.87e9,
    mm39 = 1.87e9
  )
  size <- table[[tolower(genome)]]
  if (is.null(size)) {
    stop("unsupported --genome '", genome, "'. Known: ",
         paste(names(table), collapse = ", "),
         ". Add the effective genome size for your build to call_peaks.R.")
  }
  size
}

# ---------------------------------------------------------------------------
# Sample fragment barcodes once at startup (expensive on 100GB files; we
# only need a few thousand to detect the suffix pattern).
# ---------------------------------------------------------------------------

sample_fragment_barcodes <- function(path, max_rows) {
  message(sprintf("[call_peaks] sampling unique barcodes from first %s fragment rows ...",
                  format(max_rows, big.mark = ",")))
  cmd <- sprintf("gunzip -c %s | head -%d | awk '!seen[$4]++ {print $4}'",
                 shQuote(path), max_rows)
  bcs <- tryCatch(system(cmd, intern = TRUE),
                  error = function(e) {
                    stop("failed to sample fragments file via gunzip|awk: ",
                         conditionMessage(e),
                         "\nIs the file bgzipped and is gunzip on PATH?")
                  })
  message("[call_peaks]   collected ", length(bcs),
          " unique barcodes from fragments sample")
  bcs
}

# ---------------------------------------------------------------------------
# Barcode reconciliation: try a few common suffix-stripping patterns and pick
# the one with highest overlap against the fragment-barcode set. Fail loud
# if even the best correction can't reach --min-overlap.
# ---------------------------------------------------------------------------

reconcile_barcodes <- function(group_cells, frag_bcs, min_overlap, group_name) {
  raw <- length(intersect(group_cells, frag_bcs)) / length(group_cells)

  candidates <- list(
    list(name = "none",                       fn = function(x) x),
    list(name = "split _ take first",         fn = function(x) sapply(strsplit(x, "_", fixed = TRUE), `[`, 1)),
    list(name = "split - take first",         fn = function(x) sapply(strsplit(x, "-", fixed = TRUE), `[`, 1)),
    list(name = "strip trailing -<digits>",   fn = function(x) sub("-[0-9]+$", "", x)),
    list(name = "strip trailing _<digits>",   fn = function(x) sub("_[0-9]+$", "", x))
  )

  best <- list(overlap = -1, name = NA_character_, new_cells = group_cells)
  for (c in candidates) {
    new_cells <- c$fn(group_cells)
    overlap <- length(intersect(new_cells, frag_bcs)) / length(new_cells)
    if (overlap > best$overlap) {
      best <- list(overlap = overlap, name = c$name, new_cells = new_cells)
    }
  }

  message(sprintf(
    "[call_peaks]   group %s: raw overlap %.1f%%; best %.1f%% via [%s]",
    group_name, 100 * raw, 100 * best$overlap, best$name))

  if (best$overlap < min_overlap) {
    sample_obj <- head(group_cells, 5)
    sample_frag <- head(frag_bcs, 5)
    stop(sprintf(
      paste0("group %s: barcode overlap %.1f%% (best after suffix correction) ",
             "below --min-overlap %.1f%%.\n",
             "  Sample object barcodes: %s\n",
             "  Sample fragment barcodes: %s\n",
             "Likely causes: wrong fragments file for this sample, ",
             "non-suffix barcode hashing, or cells subsetted to a different ",
             "10x run. Manually align before retrying."),
      group_name, 100 * best$overlap, 100 * min_overlap,
      paste(sample_obj, collapse = ", "),
      paste(sample_frag, collapse = ", ")
    ))
  }

  list(new_cells = best$new_cells, correction = best$name, overlap = best$overlap)
}

# ---------------------------------------------------------------------------
# narrowPeak conversion (ENCODE 10-col)
# ---------------------------------------------------------------------------
#
# Signac CallPeaks returns a GRanges with: name, score, fold_change,
# neg_log10pvalue_summit, neg_log10qvalue_summit, relative_summit_position.
# narrowPeak (BED6+4) expects:
#   chrom chromStart chromEnd name score strand signalValue pValue qValue peak

peaks_to_narrowpeak_df <- function(peaks) {
  if (length(peaks) == 0) {
    return(NULL)
  }
  data.frame(
    chrom       = as.character(seqnames(peaks)),
    chromStart  = pmax(0L, start(peaks) - 1L),    # BED is 0-based half-open
    chromEnd    = end(peaks),
    name        = if (!is.null(peaks$name)) peaks$name else paste0("peak_", seq_along(peaks)),
    score       = if (!is.null(peaks$score)) as.integer(round(peaks$score)) else 0L,
    strand      = ".",
    signalValue = if (!is.null(peaks$fold_change)) peaks$fold_change else 0,
    pValue      = if (!is.null(peaks$neg_log10pvalue_summit)) peaks$neg_log10pvalue_summit else -1,
    qValue      = if (!is.null(peaks$neg_log10qvalue_summit)) peaks$neg_log10qvalue_summit else -1,
    peak        = if (!is.null(peaks$relative_summit_position)) peaks$relative_summit_position else -1,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Pick an ATAC assay from a Seurat object (handles both 'ATAC' and 'peaks'
# canonical names plus any ChromatinAssay-class instance).
# ---------------------------------------------------------------------------

pick_atac_assay <- function(obj) {
  for (name in c("ATAC", "peaks")) {
    if (name %in% Assays(obj)) return(name)
  }
  hits <- Assays(obj)[sapply(Assays(obj),
                             function(a) inherits(obj[[a]], "ChromatinAssay"))]
  if (length(hits) == 0) {
    stop("no ATAC/ChromatinAssay found; cannot call peaks. ",
         "Available assays: ", paste(Assays(obj), collapse = ", "))
  }
  hits[1]
}

# ---------------------------------------------------------------------------
# Per-group worker
# ---------------------------------------------------------------------------

call_peaks_one <- function(g, obj, clusters_df, atac_assay, frag_bcs,
                           macs_path, gsize, output_dir, min_overlap) {
  out_rds  <- file.path(output_dir, sprintf("%s_peaks.rds", g$name))
  out_peak <- file.path(output_dir, sprintf("%s.narrowPeak", g$name))

  # Idempotency: skip if both outputs are present and non-empty.
  if (file.exists(out_rds) && file.exists(out_peak) && file.size(out_peak) > 0) {
    message("[call_peaks] group ", g$name,
            " already has outputs; skipping.")
    return(out_peak)
  }

  # Resolve the cell set for this group from clusters.csv. Group spec
  # carries metadata as a named list/dict {col: val, ...} (or NULL when
  # no metadata stratification was requested); compose the filter from
  # all key-value pairs.
  cl <- as.character(g$cluster)
  metadata <- g$metadata
  sel <- as.character(clusters_df$seurat_cluster) == cl
  if (!is.null(metadata) && length(metadata) > 0) {
    for (mcol in names(metadata)) {
      if (!(mcol %in% colnames(clusters_df))) {
        warning("[call_peaks] group ", g$name,
                ": metadata col '", mcol, "' not in clusters.csv; skipping.")
        return(NULL)
      }
      sel <- sel &
             !is.na(clusters_df[[mcol]]) &
             as.character(clusters_df[[mcol]]) == as.character(metadata[[mcol]])
    }
  }
  cells <- clusters_df$barcode[sel]

  # Restrict to ATAC-side cells in a co-embedded object, when an `assay`
  # meta.data flag exists. Otherwise use cells present in the ATAC assay.
  obj_md <- obj@meta.data
  obj_md$.barcode <- rownames(obj_md)
  if (opt$`assay-col` %in% colnames(obj_md)) {
    atac_md <- obj_md[obj_md[[opt$`assay-col`]] %in% c("ATAC", "atac") &
                      obj_md$.barcode %in% cells, , drop = FALSE]
    cells <- atac_md$.barcode
  } else {
    cells <- intersect(cells, Cells(obj[[atac_assay]]))
  }
  if (length(cells) == 0) {
    warning("[call_peaks] group ", g$name,
            ": no ATAC cells match; skipping.")
    return(NULL)
  }

  message(sprintf("[call_peaks] group %s: %d ATAC cells",
                  g$name, length(cells)))

  # Reconcile barcodes against fragments before subsetting (so we know the
  # rename map up-front).
  rec <- reconcile_barcodes(cells, frag_bcs, min_overlap, g$name)

  # Subset → restrict to ATAC assay → rename to fragment-aligned barcodes.
  obj_sub <- subset(obj, cells = cells)
  DefaultAssay(obj_sub) <- atac_assay
  # Rename only the subset object (avoids duplicate-name collisions when
  # both RNA and ATAC sides of a co-embedded object share bare barcodes).
  obj_sub <- RenameCells(obj_sub, new.names = rec$new_cells)

  # Attach a fragment object scoped to just this group's cells.
  frags <- CreateFragmentObject(path = opt$fragments, cells = rec$new_cells)
  Fragments(obj_sub) <- NULL
  Fragments(obj_sub) <- frags

  # Tag everyone in the subset with one ident so Signac::CallPeaks emits a
  # single peak set for this group.
  Idents(obj_sub) <- factor(rep(g$name, ncol(obj_sub)))

  message("[call_peaks] group ", g$name, ": running Signac::CallPeaks ...")
  peaks <- tryCatch(
    Signac::CallPeaks(
      object                = obj_sub,
      macs2.path            = macs_path,
      effective.genome.size = gsize,
      combine.peaks         = TRUE
    ),
    error = function(e) {
      warning("[call_peaks] group ", g$name, ": CallPeaks failed: ",
              conditionMessage(e), " — leaving narrowPeak unwritten.")
      NULL
    }
  )

  if (is.null(peaks) || length(peaks) == 0) {
    warning("[call_peaks] group ", g$name, ": no peaks returned.")
    return(NULL)
  }

  saveRDS(peaks, out_rds)
  df <- peaks_to_narrowpeak_df(peaks)
  if (is.null(df)) {
    warning("[call_peaks] group ", g$name, ": empty peak set; not writing narrowPeak.")
    return(NULL)
  }
  write.table(df, file = out_peak, sep = "\t",
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  message("[call_peaks] group ", g$name, ": wrote ", nrow(df),
          " peaks -> ", basename(out_peak))
  out_peak
}

# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

message("[call_peaks] loading object ", opt$input)
obj <- load_input(opt$input)
atac_assay <- pick_atac_assay(obj)
message("[call_peaks] using ATAC assay: ", atac_assay)

clusters_df <- read.csv(opt$clusters, stringsAsFactors = FALSE,
                        check.names = FALSE)
stopifnot("barcode" %in% colnames(clusters_df),
          "seurat_cluster" %in% colnames(clusters_df))

groups_plan <- fromJSON(opt$groups, simplifyDataFrame = FALSE)

frag_bcs <- sample_fragment_barcodes(opt$fragments, opt$`frag-sample`)
gsize <- effective_genome_size(opt$genome)

written <- list()
for (g in groups_plan$groups) {
  out <- tryCatch(
    call_peaks_one(g, obj, clusters_df, atac_assay, frag_bcs,
                   opt$`macs-path`, gsize, opt$`output-dir`, opt$`min-overlap`),
    error = function(e) {
      warning("[call_peaks] group ", g$name, ": ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(out)) {
    written[[length(written) + 1]] <- list(name = g$name, path = out)
  }
}

writeLines(
  toJSON(written, pretty = TRUE, auto_unbox = TRUE),
  file.path(opt$`output-dir`, "_index.json")
)

message("[call_peaks] done. ", length(written), " / ",
        length(groups_plan$groups), " groups produced narrowPeaks.")
