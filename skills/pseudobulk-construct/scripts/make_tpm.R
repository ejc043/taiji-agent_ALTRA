#!/usr/bin/env Rscript
# make_tpm.R
#
# Convert raw-count RNA pseudobulk TSVs (aggregate_rna.R output) to TPM.
# Promoted from runs/MARCH/code/make_tpm.R + runs/cima_full/code/make_tpm.R.
#
# Gene length = union-exon length per gene (exons of all transcripts merged
# with GenomicRanges::reduce, widths summed), taken from the GTF.
#
# TPM = (count / length) / sum(count / length) * 1e6, computed over the genes
# with a known length. This is the same formula as
# DGEobj.utils::convertCounts(unit = "TPM"), inlined to avoid that dependency.
#
# Gene key: aggregate_rna.R writes rownames(counts), which are usually gene
# symbols but can be Ensembl IDs. Lengths are built for both gene_name and
# gene_id (version suffix stripped), and the key matching more of the input
# genes is used.
#
# Unmatched genes (no length in the GTF) are written with TPM 0 and counted in
# tpm_summary.json. Fails if fewer than --min-match-frac of genes match, which
# almost always means the GTF is the wrong genome/annotation.
#
# Reads:
#   --gtf          GTF (plain or .gz)
#   --counts-dir   directory of <group>.tsv (gene<TAB>count, no header)
#   --output-dir   where to write <group>.tsv (gene<TAB>tpm, no header)
#
# Writes:
#   <output-dir>/<group>.tsv
#   <output-dir>/tpm_summary.json

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(rtracklayer)
  library(jsonlite)
})

option_list <- list(
  make_option("--gtf", type = "character"),
  make_option("--counts-dir", type = "character"),
  make_option("--output-dir", type = "character"),
  make_option("--min-match-frac", type = "double", default = 0.5,
              help = "Fail if fewer than this fraction of genes match the GTF. Default: 0.5.")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (k in c("gtf", "counts-dir", "output-dir")) {
  if (is.null(opt[[k]])) stop("--", k, " is required")
}
if (!file.exists(opt$gtf)) stop("GTF not found: ", opt$gtf)
dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)

message("[make_tpm] loading GTF: ", opt$gtf)
gtf <- rtracklayer::import(opt$gtf)
exons <- gtf[gtf$type == "exon"]
rm(gtf); invisible(gc())
message("[make_tpm]   exons: ", length(exons))

union_exon_lengths <- function(key) {
  keep <- !is.na(key) & key != ""
  by_gene <- split(exons[keep], key[keep])
  lens <- sum(width(reduce(by_gene)))
  lens[lens > 0]
}
len_by_name <- if ("gene_name" %in% colnames(mcols(exons)))
  union_exon_lengths(exons$gene_name) else integer(0)
len_by_id <- if ("gene_id" %in% colnames(mcols(exons)))
  union_exon_lengths(sub("\\.[0-9]+$", "", exons$gene_id)) else integer(0)
message("[make_tpm]   lengths: ", length(len_by_name), " gene_name, ",
        length(len_by_id), " gene_id")

tsvs <- list.files(opt$`counts-dir`, pattern = "\\.tsv$", full.names = TRUE)
if (length(tsvs) == 0) stop("no .tsv files in ", opt$`counts-dir`)
message("[make_tpm] converting ", length(tsvs), " pseudobulk TSVs")

# Pick the gene key once, from the first file (all files share rownames).
first_genes <- read.delim(tsvs[1], header = FALSE, stringsAsFactors = FALSE)[[1]]
frac_name <- mean(first_genes %in% names(len_by_name))
frac_id   <- mean(sub("\\.[0-9]+$", "", first_genes) %in% names(len_by_id))
use_id <- frac_id > frac_name
gene_lengths <- if (use_id) len_by_id else len_by_name
key_fn <- if (use_id) function(g) sub("\\.[0-9]+$", "", g) else identity
match_frac <- max(frac_name, frac_id)
message(sprintf("[make_tpm]   matching on %s: %.1f%% of %d genes have a length",
                if (use_id) "gene_id" else "gene_name",
                100 * match_frac, length(first_genes)))
if (match_frac < opt$`min-match-frac`) {
  stop(sprintf(paste0("only %.1f%% of genes matched the GTF (min %.0f%%). ",
                      "Is --gtf the same genome/annotation used for the ",
                      "single-cell counts?"),
               100 * match_frac, 100 * opt$`min-match-frac`))
}

per_file <- list()
for (f in tsvs) {
  df <- read.delim(f, header = FALSE, col.names = c("gene", "count"),
                   stringsAsFactors = FALSE)
  lens <- unname(gene_lengths[key_fn(df$gene)])
  matched <- !is.na(lens)
  rate <- ifelse(matched, df$count / lens, 0)
  total <- sum(rate)
  tpm <- if (total > 0) rate / total * 1e6 else rate
  write.table(data.frame(gene = df$gene, tpm = round(tpm, 4)),
              file.path(opt$`output-dir`, basename(f)), sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  per_file[[basename(f)]] <- list(
    n_genes = nrow(df),
    n_unmatched = sum(!matched),
    unmatched_count_frac = if (sum(df$count) > 0)
      round(sum(df$count[!matched]) / sum(df$count), 4) else 0
  )
  if (length(per_file) %% 100 == 0) {
    message("[make_tpm]   ", length(per_file), " / ", length(tsvs))
  }
}

writeLines(
  toJSON(list(
    gtf = normalizePath(opt$gtf),
    gene_key = if (use_id) "gene_id" else "gene_name",
    gene_length = "union-exon",
    match_frac = round(match_frac, 4),
    n_files = length(per_file),
    files = per_file
  ), pretty = TRUE, auto_unbox = TRUE),
  file.path(opt$`output-dir`, "tpm_summary.json")
)
message("[make_tpm] wrote ", length(per_file), " TPM TSVs to ", opt$`output-dir`)
