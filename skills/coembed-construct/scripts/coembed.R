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
  # GenomeInfoDb provides `seqlevelsStyle<-` (the assignment generic used to
  # flip Ensembl-style "1, 2, 3" chromosome names to UCSC-style "chr1, chr2,
  # chr3" on the EnsDb annotations). GenomicRanges DEPENDS on it but in some
  # bioconda builds the assignment form is not re-exported under
  # suppressPackageStartupMessages, so we attach it explicitly.
  library(GenomeInfoDb)
  library(jsonlite)
})

# macOS R imposes a default VSIZE cap (~16-18 GB on Apple Silicon, regardless
# of physical RAM). FindTransferAnchors / CCA / joint PCA on 60k+21k cells
# easily blow past it, errors out with "vector memory limit of 18.0 Gb
# reached, see mem.maxVSize()". On Linux there is no such cap, so this
# is a no-op there. We bump to 128 GB unconditionally — if the host has
# less physical RAM, the kernel OOMs as usual; we just want R itself to
# stop self-throttling. Override via TAIJI_R_MAX_VSIZE_GB env var.
if (Sys.info()[["sysname"]] == "Darwin") {
  cap_gb <- suppressWarnings(as.numeric(Sys.getenv("TAIJI_R_MAX_VSIZE_GB",
                                                    unset = "128")))
  if (is.na(cap_gb) || cap_gb <= 0) cap_gb <- 128
  tryCatch({
    mem.maxVSize(cap_gb * 1024)   # mem.maxVSize takes MB
    message(sprintf(
      "[coembed] macOS detected; bumped R VSIZE cap to %.0f GB (override via TAIJI_R_MAX_VSIZE_GB).",
      cap_gb))
  }, error = function(e) {
    message("[coembed] WARN: could not bump R VSIZE cap: ",
            conditionMessage(e))
  })
}

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
  make_option("--reuse-rna-reductions", action = "store_true", default = FALSE,
              help = "If RNA object already has 'pca' + 'umap' reductions and 'data' slot is populated, skip NormalizeData/FindVariableFeatures/ScaleData/RunPCA/RunUMAP. Saves ~5-10 min on 60k+ cells."),
  make_option("--reuse-atac-reductions", action = "store_true", default = FALSE,
              help = "If ATAC object already has 'lsi' + 'umap.atac' (or 'umap') reductions, skip TF-IDF/SVD/UMAP. Saves ~3-5 min on 20k+ cells."),
  make_option("--strict-metadata", action = "store_true", default = FALSE,
              help = "Fail-loud if any --metadata-cols values differ between the RNA and ATAC inputs (e.g. tissue=spleen on RNA, tissue=Spleen on ATAC creates 2 separate groups downstream instead of 1). Default: warn but continue."),
  make_option("--fragments",           type = "character", default = NULL,
              help = "Path to fragments.tsv.gz to attach to the ATAC ChromatinAssay. Required when the path embedded in the input .rds is stale (very common for objects shipped from HPC: the embedded path points at the original user's home dir and file.exists() returns FALSE on the destination machine). The skill rebuilds the Fragments handle in place using the cells already attached to the existing fragment objects, so no barcode reconciliation happens here. If the embedded path resolves on the current machine, this flag is unused."),
  make_option("--no-cluster",          action = "store_true", default = FALSE),
  make_option("--no-plot",             action = "store_true", default = FALSE),
  make_option("--abort-on-memory-risk", action = "store_true", default = FALSE,
              help = "Refuse to run if the macOS pre-flight estimates peak RAM > available budget. Default: warn but continue. Use this on shared Macs where you'd rather fail-fast than risk an OOM-kill mid-pipeline.")
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
# Repair stale fragment paths
# ------------------------------------------------------------------------
# A ChromatinAssay carries Fragment objects that hard-code an absolute
# path to fragments.tsv.gz. When the .rds is built on HPC and shipped to
# a user's laptop (very common), that path points at the original user's
# home dir and resolves nowhere — but readRDS succeeds, so the failure
# only surfaces ~10 min later inside GeneActivity() with an opaque
# tabix/seek error. Catch this at startup: if the embedded path is
# missing AND a --fragments override was provided, rebuild a single
# Fragment object in place using the existing barcode-cell mapping; if
# the path is missing AND no override, refuse loudly with the fix.
repair_fragments <- function(atac, atac_assay, fragments_override = NULL) {
  fr <- Fragments(atac[[atac_assay]])
  if (length(fr) == 0) {
    stop("ATAC ChromatinAssay '", atac_assay, "' has no Fragment objects. ",
         "GeneActivity() requires a fragments handle. Attach one upstream ",
         "via Fragments(atac[['", atac_assay, "']]) <- CreateFragmentObject(...).")
  }
  embedded_paths <- vapply(fr, function(f) tryCatch(slot(f, "path"),
                                                     error = function(e) ""),
                            character(1))
  exists_mask <- file.exists(embedded_paths)
  message(sprintf("[coembed]   ATAC fragments: %d handle(s); paths exist=%s",
                  length(fr), paste(exists_mask, collapse = ",")))
  for (i in seq_along(fr)) {
    message(sprintf("[coembed]     [%d] %s -> %s",
                    i, embedded_paths[i],
                    if (exists_mask[i]) "OK" else "MISSING"))
  }
  if (all(exists_mask)) return(atac)

  if (is.null(fragments_override)) {
    stop("ATAC fragments path embedded in the .rds does not resolve on this ",
         "machine: ", paste(embedded_paths[!exists_mask], collapse = ", "),
         ". This is the typical 'object built on HPC, scp'd to laptop' wart. ",
         "Pass --fragments <path/to/fragments.tsv.gz> so the skill can ",
         "rebuild the ChromatinAssay's fragment handle in place. ",
         "(Or fix upstream: ",
         "fr <- CreateFragmentObject(path, cells = colnames(atac)); ",
         "Fragments(atac[['", atac_assay, "']]) <- NULL; ",
         "Fragments(atac[['", atac_assay, "']]) <- fr; saveRDS(...).)")
  }

  if (!file.exists(fragments_override)) {
    stop("--fragments path '", fragments_override, "' does not exist.")
  }
  tbi <- paste0(fragments_override, ".tbi")
  if (!file.exists(tbi)) {
    message("[coembed]   WARN: tabix index '", tbi, "' not found next to ",
            "the fragments file. GeneActivity() needs it; if missing, run ",
            "`tabix -p bed ", fragments_override, "` upstream.")
  }
  # Reuse the existing cell list. If multiple Fragment objects existed
  # (rare for separate-assay), unify them to a single override path with
  # the union of cells, so the rebuilt handle is unambiguous.
  cells <- unique(unlist(lapply(fr, function(f) tryCatch(slot(f, "cells"),
                                                          error = function(e) NULL))))
  if (length(cells) == 0) {
    cells <- colnames(atac[[atac_assay]])
  }
  message(sprintf("[coembed]   rebuilding Fragment handle: path=%s, cells=%d",
                  fragments_override, length(cells)))
  new_fr <- CreateFragmentObject(path = fragments_override,
                                 cells = cells, validate.fragments = FALSE,
                                 verbose = FALSE)
  Fragments(atac[[atac_assay]]) <- NULL
  Fragments(atac[[atac_assay]]) <- new_fr
  atac
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

validate_metadata_consistency <- function(rna, atac, meta_cols, strict = FALSE) {
  # Catches the silent-bug case where the same conceptual metadata col has
  # different value SETS or different value CASING across the two inputs.
  #
  # Real example: RNA has tissue ∈ {spleen, siIEL}, ATAC has tissue ∈
  # {Spleen, siIEL}. After merge the column has THREE distinct values
  # (siIEL appears once in each, but spleen + Spleen are split). Then
  # pseudobulk-construct --metadata-cols tissue produces 3 cross-product
  # groups per cluster instead of 2 — silently doubling the spleen
  # samples and giving each half-strength. Loudly warn here so the user
  # harmonizes upstream before spending 30 min on coembed.
  rna_md  <- attributes(rna)$meta.data
  atac_md <- attributes(atac)$meta.data
  any_issue <- FALSE
  for (col in meta_cols) {
    on_rna  <- col %in% colnames(rna_md)
    on_atac <- col %in% colnames(atac_md)
    if (!on_rna && !on_atac) {
      message(sprintf("[coembed] WARN: --metadata-cols '%s' is on neither input", col))
      any_issue <- TRUE
      next
    }
    if (!on_rna || !on_atac) {
      message(sprintf("[coembed] WARN: --metadata-cols '%s' on only %s; cells from the other side will be NA",
                      col, if (on_rna) "RNA" else "ATAC"))
      any_issue <- TRUE
      next
    }
    rna_vals  <- sort(unique(as.character(rna_md[[col]])))
    atac_vals <- sort(unique(as.character(atac_md[[col]])))
    rna_only  <- setdiff(rna_vals, atac_vals)
    atac_only <- setdiff(atac_vals, rna_vals)
    if (length(rna_only) > 0 || length(atac_only) > 0) {
      # Check for case-only mismatch — most common silent bug.
      ci_match <- length(setdiff(tolower(rna_vals), tolower(atac_vals))) == 0 &&
                  length(setdiff(tolower(atac_vals), tolower(rna_vals))) == 0
      severity <- if (ci_match) "CASE-ONLY" else "VALUE"
      message(sprintf(
        "[coembed] WARN: %s mismatch in --metadata-cols '%s'.\n  RNA values:  {%s}\n  ATAC values: {%s}\n  RNA-only:    {%s}\n  ATAC-only:   {%s}\n  Effect: cross-product stratification downstream will create separate groups for each variant. Harmonize upstream (e.g. atac$%s <- tolower(atac$%s)) before running coembed if you want them merged.",
        severity, col,
        paste(rna_vals, collapse=", "), paste(atac_vals, collapse=", "),
        paste(rna_only, collapse=", "), paste(atac_only, collapse=", "),
        col, col))
      any_issue <- TRUE
    }
  }
  if (any_issue && strict) {
    stop("[coembed] metadata-cols validation failed under strict mode. ",
         "Either harmonize the values upstream or remove --strict-metadata.")
  }
  invisible(any_issue)
}

binary_search_resolution <- function(obj, n_cells, target_size, max_iter = 8L) {
  r <- seed_resolution(n_cells, target_size)
  lo <- 0.01; hi <- 5.0
  trace <- list()
  obj_try <- obj
  prev_n_clust <- NA_integer_
  prev_r <- NA_real_
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
    # Saturation early-stop: if Louvain hits its modularity ceiling at the
    # outer bracket (r == hi == 5.0) and the cluster count is unchanged from
    # the prior iter, further iters at the same resolution are guaranteed to
    # produce the same partition (FindClusters is deterministic at fixed
    # seed). Break out instead of burning ~60s/iter on identical work.
    # Common when target_size is small relative to the data's intrinsic
    # cluster structure (e.g. target=200 on 81k cells where Louvain
    # plateaus at ~75 clusters / 1080 cells/cluster).
    if (!is.na(prev_n_clust) && n_clust == prev_n_clust && r == prev_r &&
        (r >= hi - 1e-9 || r <= lo + 1e-9)) {
      message(sprintf(
        "[coembed] resolution search hit Louvain modularity ceiling at r=%.3f (n_clust=%d unchanged across two iters at the bracket bound). Stopping early; downstream stratification still works on the resulting clusters.",
        r, n_clust))
      return(list(obj = obj_try, resolution = r, trace = trace))
    }
    prev_n_clust <- n_clust
    prev_r <- r
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

# ------------------------------------------------------------------------
# Seurat version compatibility check
# ------------------------------------------------------------------------
# The skill targets Seurat v5+. Inputs saved with Assay5 (Seurat v5) load
# without error in Seurat v4 BUT silently misbehave: Assays() returns
# empty, GetAssayData() errors with "'RNA' is not an assay", and the
# whole pipeline crashes minutes in on the resulting NULL counts. Catch
# this at startup so the user fixes their env before spending 30 min on
# coembed.
seurat_version <- tryCatch(
  utils::packageVersion("Seurat"),
  error = function(e) NULL
)
if (is.null(seurat_version)) {
  stop("Seurat is not installed. Run `bin/install.sh --profile sc` first.")
}
message(sprintf("[coembed] Seurat version: %s", as.character(seurat_version)))
if (seurat_version < "5.0.0") {
  message(sprintf(
    "[coembed] WARN: Seurat %s detected (< 5.0). If your input objects use ",
    as.character(seurat_version)),
    "Seurat v5 layouts (Assay5 with layered counts/data) the skill will ",
    "fail because v4 can't access Assay5 contents. Consider upgrading ",
    "the env: `micromamba install -n taiji-agent -c bioconda 'r-seurat>=5.0'`. ",
    "Continuing — if your inputs are v3 Assay objects this will work fine.")
}

message("[coembed] loading RNA from ", opt$rna)
rna <- readRDS(opt$rna)
if (!inherits(rna, "Seurat")) stop("RNA file is not a Seurat object")

# Detect Assay5 inputs under Seurat <5 and refuse loudly. If the user
# squelched the version warning above, this is the explicit guard at the
# point of failure. We check via the assay class, not Assays(rna),
# because Assays() returns empty when v4 reads a v5 object.
for (an in names(attributes(rna)$assays)) {
  ac <- class(attributes(rna)$assays[[an]])
  if ("Assay5" %in% ac && seurat_version < "5.0.0") {
    stop(sprintf(
      "RNA assay '%s' is class Assay5 (Seurat v5) but Seurat %s is loaded. ",
      an, as.character(seurat_version)),
      "Assays() returns empty under v4 and the pipeline crashes ~10 min in. ",
      "Either: (1) upgrade the env to Seurat >= 5, or (2) downgrade the ",
      "input via `obj[['RNA']] <- as(obj[['RNA']], 'Assay'); saveRDS(...)`.")
  }
}
DefaultAssay(rna) <- "RNA"

# NOTE: SeuratObject >= 5.0 made the `slot=` arg of GetAssayData() defunct
# (not just deprecated — it errors out). Use `layer=` everywhere.
rna_counts <- tryCatch(GetAssayData(rna, layer = "counts"),
                       error = function(e) NULL)
if (is.null(rna_counts) || nrow(rna_counts) == 0) {
  stop("RNA object has no raw counts. Cannot proceed with NormalizeData.")
}

# Diagnostic snapshot of RNA input — surface what's already there so the
# user knows what's about to be re-computed (or reused via --reuse-*).
message(sprintf(
  "[coembed]   RNA input: %d cells, assays={%s}, reductions={%s}, has_seurat_clusters=%s",
  ncol(rna),
  paste(Assays(rna), collapse = ","),
  paste(names(rna@reductions), collapse = ","),
  "seurat_clusters" %in% colnames(rna@meta.data)
))

# Preserve the user's per-modality cluster labels under a non-clobbering
# name BEFORE we run FindClusters on the joint PCA later. This means
# coembed@meta.data$rna_input_clusters survives the merge intact, and
# the user can compare their original RNA-only clustering to the new
# joint clustering downstream.
if ("seurat_clusters" %in% colnames(rna@meta.data)) {
  rna$rna_input_clusters <- as.character(rna$seurat_clusters)
  message("[coembed]   preserving RNA's existing seurat_clusters as 'rna_input_clusters'")
}

# Reuse path: skip the expensive standard preprocessing if the user
# asserts the input is already prepared and the slots/reductions exist.
rna_has_pca  <- "pca"  %in% names(rna@reductions)
rna_has_umap <- "umap" %in% names(rna@reductions)
rna_data_layer <- tryCatch(GetAssayData(rna, layer = "data"),
                           error = function(e) NULL)
rna_has_data <- !is.null(rna_data_layer) && nrow(rna_data_layer) > 0
rna_has_vf   <- length(VariableFeatures(rna)) > 0

if (opt$`reuse-rna-reductions` && rna_has_pca && rna_has_umap &&
    rna_has_data && rna_has_vf) {
  message("[coembed]   --reuse-rna-reductions: keeping existing pca/umap/data/VFs.")
} else {
  if (opt$`reuse-rna-reductions`) {
    message("[coembed]   --reuse-rna-reductions set but RNA object missing prerequisites ",
            "(pca=", rna_has_pca, ", umap=", rna_has_umap,
            ", data slot=", rna_has_data, ", variable features=", rna_has_vf,
            "). Falling back to full preprocessing.")
  }
  message("[coembed] RNA: standard preprocessing (Norm/VF/Scale/PCA/UMAP)")
  rna <- NormalizeData(rna, verbose = FALSE)
  rna <- FindVariableFeatures(rna, verbose = FALSE)
  rna <- ScaleData(rna, verbose = FALSE)
  rna <- RunPCA(rna, npcs = opt$`n-pcs`, verbose = FALSE)
  rna <- RunUMAP(rna, dims = 1:opt$`n-pcs`, verbose = FALSE)
}
rna$assay <- "RNA"

# ------------------------------------------------------------------------
# Stage 2 — load & process ATAC, set annotation, RunTFIDF/SVD/UMAP
# ------------------------------------------------------------------------

message("[coembed] loading ATAC from ", opt$atac)
atac <- readRDS(opt$atac)
if (!inherits(atac, "Seurat")) stop("ATAC file is not a Seurat object")

# Same Assay5-under-v4 guard for ATAC.
for (an in names(attributes(atac)$assays)) {
  ac <- class(attributes(atac)$assays[[an]])
  if ("Assay5" %in% ac && seurat_version < "5.0.0") {
    stop(sprintf(
      "ATAC assay '%s' is class Assay5 (Seurat v5) but Seurat %s is loaded. ",
      an, as.character(seurat_version)),
      "Upgrade the env to Seurat >= 5, or downgrade the input.")
  }
}

atac_assay <- pick_atac_assay(atac)
DefaultAssay(atac) <- atac_assay
message("[coembed] ATAC: using assay '", atac_assay, "'")

# Diagnostic snapshot.
message(sprintf(
  "[coembed]   ATAC input: %d cells, assays={%s}, reductions={%s}, has_seurat_clusters=%s",
  ncol(atac),
  paste(Assays(atac), collapse = ","),
  paste(names(atac@reductions), collapse = ","),
  "seurat_clusters" %in% colnames(atac@meta.data)
))

# Repair stale fragments path BEFORE GeneActivity / annotation. If the
# .rds was built on a different machine, the embedded path in each
# Fragment object is absolute and will not resolve here. GeneActivity()
# would fail ~10 min in with an opaque tabix error.
atac <- repair_fragments(atac, atac_assay,
                         fragments_override = opt$fragments)

# Some ATAC objects (e.g. those from prior integrate_atac runs, or split
# multiome objects) already carry an "RNA" assay. We're going to overwrite
# it with freshly imputed values from TransferData later — surface this
# loudly so the user isn't surprised that imputation replaces what was
# there.
if ("RNA" %in% Assays(atac)) {
  message("[coembed]   NOTE: ATAC object already carries an 'RNA' assay. ",
          "It will be REPLACED later by TransferData-imputed values from ",
          "the RNA reference. If you want to keep the existing values, ",
          "save them under a different assay name before running this skill.")
}

# Preserve ATAC's existing seurat_clusters under a non-clobbering name.
if ("seurat_clusters" %in% colnames(atac@meta.data)) {
  atac$atac_input_clusters <- as.character(atac$seurat_clusters)
  message("[coembed]   preserving ATAC's existing seurat_clusters as 'atac_input_clusters'")
}

# ------------------------------------------------------------------------
# Cross-input metadata sanity check — catch capitalization / value-set
# mismatches BEFORE the long-running steps fire.
# ------------------------------------------------------------------------
if (length(meta_cols) > 0) {
  validate_metadata_consistency(rna, atac, meta_cols,
                                strict = opt$`strict-metadata`)
}

# ------------------------------------------------------------------------
# Memory pre-flight on macOS — catch the swap-thrash trap BEFORE GeneActivity
# ------------------------------------------------------------------------
# On 18 GB-physical Macs (which is the most common Mac config in our lab),
# FindTransferAnchors + joint PCA peaks at ~10-12 GB RSS for a 60k+21k pair.
# The default macOS swap budget is ~24 GB, and the spike fits — but barely.
# For larger pairs (e.g. 80k+30k) the working set exceeds RAM+swap and the
# kernel goes into thrash or jetsam-kills R mid-CCA with no R error message.
#
# We compute a pessimistic peak-RSS estimate from cell counts (with fudge=6
# for L2-norm / score / scaled-data intermediates) and emit:
#   - INFO when peak < 0.8 * physical_RAM (everything fits in RAM)
#   - WARN when peak ∈ [0.8 * RAM, RAM + 0.5 * swap] (slow but should finish)
#   - HARD-WARN when peak > RAM + 0.5 * swap (likely OOM-kill on macOS)
#
# Linux is unaffected (no jetsam, no fixed swap budget, R doesn't self-cap).
# We still log the estimate so SLURM users can sanity-check their --mem.
n_features_planned <- max(length(VariableFeatures(rna)), 2000L)
peak_gb_estimate <- 8 * n_features_planned * (ncol(rna) + ncol(atac)) * 6 /
                    1024 / 1024 / 1024
message(sprintf(
  "[coembed] peak-RAM estimate (CCA + joint PCA): %.1f GB at %d cells x %d features",
  peak_gb_estimate, ncol(rna) + ncol(atac), n_features_planned))

if (Sys.info()[["sysname"]] == "Darwin") {
  ram_gb <- tryCatch({
    as.numeric(system2("sysctl", c("-n", "hw.memsize"), stdout = TRUE)) /
      1024 / 1024 / 1024
  }, error = function(e) NA_real_)
  swap_free_gb <- tryCatch({
    line <- system2("sysctl", c("-n", "vm.swapusage"), stdout = TRUE)
    if (grepl("free = ", line, fixed = TRUE)) {
      as.numeric(sub(".*free = ([0-9.]+)M.*", "\\1", line)) / 1024
    } else NA_real_
  }, error = function(e) NA_real_)

  if (!is.na(ram_gb)) {
    budget_gb <- ram_gb + 0.5 * ifelse(is.na(swap_free_gb), 0, swap_free_gb)
    message(sprintf(
      "[coembed]   macOS budget: %.0f GB RAM + %.1f GB free swap = %.1f GB usable",
      ram_gb, ifelse(is.na(swap_free_gb), 0, swap_free_gb), budget_gb))
    if (peak_gb_estimate > budget_gb) {
      message(sprintf(
        "[coembed]   HARD-WARN: estimated peak %.1f GB > usable budget %.1f GB.\n  This run is very likely to be OOM-killed by macOS jetsam during\n  FindTransferAnchors (~5 min in) or joint PCA (~15 min in), with NO\n  R error reported. Recommended mitigations (any one usually unblocks):\n    1. Run on SLURM (`sbatch --mem=128G ...`) — the skill is unchanged on Linux.\n    2. Downsample the RNA reference: `rna <- subset(rna, downsample = 20000); saveRDS(rna, ...)`.\n       ATAC stays full-size (it's the query, not the working-set driver).\n    3. Lower --n-pcs from 30 to 20 (smaller joint PCA scaled matrix).\n    4. Free swap before running: close browsers / quit big GUI apps.\n  Continuing — pass --abort-on-memory-risk to refuse instead.",
        peak_gb_estimate, budget_gb))
      if (isTRUE(opt$`abort-on-memory-risk`)) {
        stop(sprintf(
          "[coembed] aborting: estimated peak %.1f GB exceeds macOS budget %.1f GB and --abort-on-memory-risk is set.",
          peak_gb_estimate, budget_gb))
      }
    } else if (peak_gb_estimate > 0.8 * ram_gb) {
      message(sprintf(
        "[coembed]   WARN: estimated peak %.1f GB > 0.8 * physical RAM (%.1f GB). Run will fit but rely on swap; expect 2-3x slower CCA.",
        peak_gb_estimate, ram_gb))
    }
  }
}

# Annotation pull for GeneActivity. Signac wants UCSC-style chromosome names
# and a `genome` attribute on the GRanges; this is what the vignette does.
ensdb <- load_ensdb(opt$genome)
message("[coembed] pulling gene annotations from EnsDb (", opt$genome, ") ...")
annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
# seqlevelsStyle<- may attempt a UCSC network fetch on HPC nodes without
# internet. Fall back to manual chr-prefix renaming if it fails.
tryCatch(
  seqlevelsStyle(annotations) <- "UCSC",
  error = function(e) {
    message("[coembed]   WARN: seqlevelsStyle UCSC network fetch failed (", conditionMessage(e), ").")
    message("[coembed]   Falling back to manual Ensembl -> UCSC chromosome rename.")
    lvls <- seqlevels(annotations)
    new_lvls <- ifelse(grepl("^chr", lvls), lvls,
                       ifelse(lvls == "MT", "chrM", paste0("chr", lvls)))
    seqlevels(annotations) <<- new_lvls
  }
)
genome(annotations) <- opt$genome
Annotation(atac) <- annotations

# Reuse path for ATAC: skip TF-IDF/SVD/UMAP if reductions already there.
# We accept either umap.atac (the Stuart vignette name) or the bare
# umap (since some users save it that way) for the existing-UMAP check.
atac_has_lsi  <- "lsi" %in% names(atac@reductions)
atac_has_umap <- any(c("umap.atac", "umap") %in% names(atac@reductions))
lsi_dims <- (opt$`lsi-skip-first` + 1):opt$`n-pcs`

if (opt$`reuse-atac-reductions` && atac_has_lsi && atac_has_umap) {
  message("[coembed]   --reuse-atac-reductions: keeping existing lsi/umap.")
  # If the umap reduction is named bare 'umap', alias it to 'umap.atac' so
  # the pre-merge plot can find it under a stable name. We don't alter the
  # underlying reduction; just add a name-pointer.
  if (!"umap.atac" %in% names(atac@reductions) &&
      "umap" %in% names(atac@reductions)) {
    atac@reductions[["umap.atac"]] <- atac@reductions[["umap"]]
  }
} else {
  if (opt$`reuse-atac-reductions`) {
    message("[coembed]   --reuse-atac-reductions set but ATAC object missing ",
            "prerequisites (lsi=", atac_has_lsi, ", umap=", atac_has_umap,
            "). Falling back to full preprocessing.")
  }
  message("[coembed] ATAC: TF-IDF / FindTopFeatures / RunSVD / RunUMAP")
  atac <- RunTFIDF(atac, verbose = FALSE)
  atac <- FindTopFeatures(atac, min.cutoff = "q0", verbose = FALSE)
  atac <- RunSVD(atac, verbose = FALSE)
  atac <- RunUMAP(atac, reduction = "lsi", dims = lsi_dims,
                  reduction.name = "umap.atac", verbose = FALSE)
}
atac$assay <- "ATAC"

# ------------------------------------------------------------------------
# Stage 3 — Gene activity assay for anchor transfer
# ------------------------------------------------------------------------
# If the ATAC prep object already has an "RNA" assay (pre-built gene
# activity, common in _prep objects from this lab's pipeline), reuse it
# directly — this matches coembed_preprocess.R which used query.assay="RNA".
# Computing fresh GeneActivity restricted to 2000 variable features drops
# ~18% of features before anchoring even starts (1634/2000 shared in D7).
# Using the full pre-built assay (21808 genes) gives FindTransferAnchors
# access to the complete variable feature set.
#
# When no pre-built RNA assay exists, compute GeneActivity as a fallback
# (stored as "ACTIVITY").

atac_gene_activity_assay <- if ("RNA" %in% Assays(atac)) {
  message("[coembed] ATAC object has pre-built 'RNA' assay (gene activity). ",
          "Reusing it for anchor transfer (matches lab reference coembed_preprocess.R).")
  DefaultAssay(atac) <- "RNA"
  atac <- NormalizeData(atac, verbose = FALSE)
  atac <- ScaleData(atac, features = rownames(atac[["RNA"]]), verbose = FALSE)
  "RNA"
} else {
  message("[coembed] computing per-cell gene activities (ATAC -> ACTIVITY assay)")
  gene_activities <- GeneActivity(atac, features = VariableFeatures(rna))
  atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene_activities)
  DefaultAssay(atac) <- "ACTIVITY"
  atac <- NormalizeData(atac, verbose = FALSE)
  atac <- ScaleData(atac, features = rownames(atac), verbose = FALSE)
  "ACTIVITY"
}

# ------------------------------------------------------------------------
# Stage 4 — FindTransferAnchors (RNA ref, ATAC query, CCA)
# ------------------------------------------------------------------------

# Features-used audit: FindTransferAnchors silently drops any feature in
# `features = ...` that isn't present in BOTH the reference and query
# assays, and only prints a one-line warning to stderr. Compute the
# intersection up-front so we can record it in the summary JSON. A high
# drop rate (>50%) usually means the genome / EnsDb / RNA gene-name
# convention don't match (e.g. running mm10 EnsDb against a human RNA
# object, or Ensembl IDs vs gene symbols).
features_requested <- VariableFeatures(rna)
features_in_atac   <- intersect(features_requested,
                                 rownames(GetAssayData(atac, assay = atac_gene_activity_assay,
                                                       layer = "data")))
n_features_used   <- length(features_in_atac)
n_features_dropped <- length(features_requested) - n_features_used
message(sprintf(
  "[coembed]   features for anchoring: %d requested -> %d shared (%.0f%% kept, %d dropped from %s assay)",
  length(features_requested), n_features_used,
  100 * n_features_used / max(1, length(features_requested)),
  n_features_dropped, atac_gene_activity_assay))
if (n_features_used < 0.5 * length(features_requested)) {
  message(sprintf(
    "[coembed]   WARN: only %d of %d RNA variable features are present in the gene activity assay. Common cause: --genome / EnsDb mismatch with the RNA gene-name convention (e.g. mm10 EnsDb on human data, or Ensembl IDs vs symbols). The Stuart vignette assumes both inputs use gene symbols matching the EnsDb's gene_name slot.",
    n_features_used, length(features_requested)))
}

message("[coembed] FindTransferAnchors: RNA (ref) -> ATAC (query, ", atac_gene_activity_assay, " assay), CCA")
transfer_anchors <- FindTransferAnchors(
  reference       = rna,
  query           = atac,
  features        = VariableFeatures(rna),
  reference.assay = "RNA",
  query.assay     = atac_gene_activity_assay,
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
refdata <- GetAssayData(rna, assay = "RNA", layer = "data")[genes_use, ]
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
# JoinLayers required in Seurat v5: merge() splits RNA into counts.1/data.1 etc;
# ScaleData/RunPCA need a single contiguous layer.
if (packageVersion("SeuratObject") >= "5.0.0") {
  coembed <- JoinLayers(coembed)
}
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
# Use umap-learn (Python) with correlation metric to match lab reference (coembed_preprocess.R).
# R/reticulate virtualenv may point to a stale Python path; override with
# RETICULATE_PYTHON so it finds the correct conda env Python that has umap-learn.
# Falls back to uwot if umap-learn is unavailable after path correction.
local({
  override <- Sys.getenv("RETICULATE_PYTHON", unset = "")
  if (nzchar(override)) {
    message("[coembed] RETICULATE_PYTHON already set: ", override)
  } else {
    # Prefer the taiji-agent conda env Python (umap-learn installed there).
    candidates <- c(
      "/stg3/data1/eunice/.local/share/mamba/envs/taiji-agent/bin/python",
      file.path(Sys.getenv("CONDA_PREFIX", ""), "bin/python")
    )
    for (p in candidates) {
      if (nzchar(p) && file.exists(p)) {
        Sys.setenv(RETICULATE_PYTHON = p)
        message("[coembed] set RETICULATE_PYTHON=", p)
        break
      }
    }
  }
})
umap_method <- tryCatch({
  reticulate::import("umap")
  "umap-learn"
}, error = function(e) {
  message("[coembed] umap-learn not available (", conditionMessage(e), "); falling back to uwot")
  "uwot"
})
message("[coembed] UMAP method: ", umap_method)
coembed <- RunUMAP(coembed, dims = 1:opt$`n-pcs`,
                   umap.method = umap_method, metric = "correlation",
                   verbose = FALSE)

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
  n_features_requested = length(features_requested),
  n_features_used      = n_features_used,
  n_features_dropped   = n_features_dropped,
  peak_ram_estimate_gb = round(peak_gb_estimate, 1),
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
# QC UMAPs — always rendered unless --no-plot
# ------------------------------------------------------------------------
# Two figures get written:
#
#   qc/umap_pre_merge.png   — Stuart-vignette-style side-by-side: RNA's
#                             own UMAP next to ATAC's own UMAP. Lets you
#                             see what each modality looked like BEFORE
#                             integration. Sanity check: clusters should
#                             look reasonable per-modality independently.
#
#   qc/umap.png             — Joint UMAP on the shared (post-merge) PCA,
#                             with one panel per group-by axis:
#                               1. seurat_clusters — the de novo joint
#                                  clusters (when not --no-cluster)
#                               2. assay — RNA vs ATAC origin, distinct
#                                  qualitative colors. The integration
#                                  QC panel: well-integrated data shows
#                                  cells of both colors mixed across
#                                  clusters; poorly-integrated data
#                                  shows clusters segregating by assay.
#                               3..N — one panel per --metadata-cols
#                                  entry (genotype, tissue, etc.).
#                                  Missing cols emit a WARN.
#
# Plus qc/umap_coords.csv with barcode + 2D coords + assay + cluster
# for re-plotting outside R (matplotlib/plotly/etc.).

if (!opt$`no-plot`) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("[coembed] ggplot2 not available; skipping QC UMAPs.")
  } else {
    qc_dir <- file.path(output_dir, "qc")
    dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

    # ----- Coords CSV (cheap, useful for downstream Python plotting) ----
    emb <- Embeddings(coembed, reduction = "umap")
    umap_df <- data.frame(
      barcode = rownames(emb),
      umap_1  = emb[, 1],
      umap_2  = emb[, 2],
      assay   = coembed$assay,
      seurat_cluster = if (!opt$`no-cluster`) as.character(Idents(coembed)) else NA,
      stringsAsFactors = FALSE
    )
    # Carry user-requested metadata cols into the CSV too, so external
    # plots can color by them without reloading the .rds.
    for (mc in intersect(meta_cols, colnames(coembed@meta.data))) {
      umap_df[[mc]] <- coembed@meta.data[[mc]]
    }
    write.csv(umap_df, file.path(qc_dir, "umap_coords.csv"), row.names = FALSE)

    # ----- Validate metadata cols up-front, warn on missing -------------
    missing_cols <- setdiff(meta_cols, colnames(coembed@meta.data))
    present_cols <- intersect(meta_cols, colnames(coembed@meta.data))
    if (length(missing_cols) > 0) {
      message("[coembed] WARN: --metadata-cols not present on merged ",
              "object (skipped panels): ",
              paste(missing_cols, collapse = ", "),
              ". These columns must exist on at least one of the input ",
              "objects to survive the merge.")
    }

    # ----- Pre-merge side-by-side (Stuart vignette style) ---------------
    pre_merge_panels <- list()
    rna_p <- tryCatch(
      DimPlot(rna, reduction = "umap", group.by = "assay") +
        ggplot2::ggtitle(sprintf("RNA only (n=%d cells)", ncol(rna))) +
        ggplot2::scale_color_manual(values = c("RNA" = "#1f77b4")) +
        Seurat::NoLegend(),
      error = function(e) {
        message("[coembed] pre-merge RNA DimPlot failed: ",
                conditionMessage(e)); NULL
      }
    )
    atac_p <- tryCatch(
      DimPlot(atac, reduction = "umap.atac", group.by = "assay") +
        ggplot2::ggtitle(sprintf("ATAC only (n=%d cells)", ncol(atac))) +
        ggplot2::scale_color_manual(values = c("ATAC" = "#ff7f0e")) +
        Seurat::NoLegend(),
      error = function(e) {
        message("[coembed] pre-merge ATAC DimPlot failed: ",
                conditionMessage(e)); NULL
      }
    )
    if (!is.null(rna_p))  pre_merge_panels[[length(pre_merge_panels) + 1]] <- rna_p
    if (!is.null(atac_p)) pre_merge_panels[[length(pre_merge_panels) + 1]] <- atac_p

    if (length(pre_merge_panels) > 0 &&
        requireNamespace("patchwork", quietly = TRUE)) {
      pre_png <- file.path(qc_dir, "umap_pre_merge.png")
      ggplot2::ggsave(
        pre_png,
        plot = patchwork::wrap_plots(pre_merge_panels, ncol = 2),
        width = 12, height = 5, dpi = 100, bg = "white"
      )
      message("[coembed] pre-merge UMAP -> ", pre_png)
    }

    # ----- Joint UMAP panels (post-merge, shared space) -----------------
    # Order: clusters (if any) -> assay (always second so it's prominent)
    #        -> each user-requested metadata col, in order they were given.
    groupby_cols <- c()
    if (!opt$`no-cluster`) groupby_cols <- c(groupby_cols, "seurat_clusters")
    groupby_cols <- c(groupby_cols, "assay", present_cols)
    groupby_cols <- unique(groupby_cols)

    # Cell counts for the title of the assay panel.
    n_rna  <- sum(coembed$assay == "RNA")
    n_atac <- sum(coembed$assay == "ATAC")

    panels <- list()
    for (col in groupby_cols) {
      title <- if (col == "assay") {
        sprintf("assay  (RNA: %d  /  ATAC: %d)", n_rna, n_atac)
      } else if (col == "seurat_clusters" && !is.na(chosen_resolution)) {
        sprintf("seurat_clusters  (res=%.3f, n=%d)",
                chosen_resolution, length(unique(Idents(coembed))))
      } else col
      p <- tryCatch({
        plot <- DimPlot(coembed, reduction = "umap", group.by = col,
                        label = (col == "seurat_clusters"), repel = TRUE) +
               ggplot2::ggtitle(title)
        # Explicit qualitative colors for the assay panel so RNA vs ATAC
        # is unmistakable on grayscale prints / colorblind viewers.
        if (col == "assay") {
          plot <- plot + ggplot2::scale_color_manual(
            values = c("RNA" = "#1f77b4", "ATAC" = "#ff7f0e"))
        }
        plot
      }, error = function(e) {
        message("[coembed] DimPlot[", col, "] failed: ", conditionMessage(e))
        NULL
      })
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
      message("[coembed] joint UMAP -> ", out_png,
              " (", length(panels), " panels: ",
              paste(groupby_cols, collapse = ", "), ")")
    }
  }
}

message("[coembed] done.")
