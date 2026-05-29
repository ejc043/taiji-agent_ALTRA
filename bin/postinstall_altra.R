#!/usr/bin/env Rscript
# bin/postinstall_altra.R
# Installs ArchR from GitHub — the only ALTRA dependency not on conda/bioconda.
# Run once after: micromamba create -f environment.altra.yml
#
# Usage:
#   micromamba activate taiji-agent-altra
#   Rscript bin/postinstall_altra.R

cat("Installing ArchR 1.0.3 from GitHub...\n")
cat("(All other dependencies were installed via environment.altra.yml)\n\n")

if (!requireNamespace("remotes", quietly = TRUE)) {
  stop("remotes package not found. Is taiji-agent-altra activated?\n",
       "Run: micromamba activate taiji-agent-altra")
}

archr_done <- requireNamespace("ArchR", quietly = TRUE)
if (archr_done) {
  cat("ArchR is already installed:", as.character(packageVersion("ArchR")), "\n")
}

# ── Compile-time fixes for source packages in the conda env ─────────────────
conda_prefix <- Sys.getenv("CONDA_PREFIX")

# 1. Point pkg-config at the conda env so harfbuzz/freetype headers are found.
if (nchar(conda_prefix) > 0) {
  conda_pc <- file.path(conda_prefix, "lib", "pkgconfig")
  current_pc <- Sys.getenv("PKG_CONFIG_PATH")
  if (!grepl(conda_pc, current_pc, fixed = TRUE)) {
    Sys.setenv(PKG_CONFIG_PATH = if (nchar(current_pc) > 0)
      paste(conda_pc, current_pc, sep = ":") else conda_pc)
    cat("PKG_CONFIG_PATH set to include", conda_pc, "\n")
  }
}

# 2. Write ~/.R/Makevars fixes needed for source compilation in the conda env:
#    a) CXX11STD: RcppArmadillo (conda) requires C++14; packages that declare
#       CXX_STD: CXX11 must be compiled with -std=gnu++14 (a strict superset).
#    b) CPPFLAGS: gcc's default search path omits $CONDA_PREFIX/include when
#       launched via `micromamba run`, so lzma.h / harfbuzz.h etc. aren't found
#       during source builds (e.g. Rhtslib bundles htslib and needs lzma.h).
makevars_dir <- file.path(Sys.getenv("HOME"), ".R")
if (!dir.exists(makevars_dir)) dir.create(makevars_dir, recursive = TRUE)
makevars_path <- file.path(makevars_dir, "Makevars")
existing <- if (file.exists(makevars_path)) readLines(makevars_path) else character(0)
flags_to_add <- list(
  CXX11STD = "CXX11STD=-std=gnu++14"
)
if (nchar(conda_prefix) > 0) {
  flags_to_add$CPPFLAGS <- paste0("CPPFLAGS += -I", conda_prefix, "/include")
  flags_to_add$LDFLAGS  <- paste0("LDFLAGS += -L", conda_prefix, "/lib")
}
for (key in names(flags_to_add)) {
  if (!any(grepl(key, existing, fixed = TRUE))) {
    write(flags_to_add[[key]], makevars_path, append = TRUE)
    cat("Added", flags_to_add[[key]], "to ~/.R/Makevars\n")
  }
}
# Also propagate include/lib paths to subprocesses (configure scripts, etc.)
if (nchar(conda_prefix) > 0) {
  inc <- file.path(conda_prefix, "include")
  lib <- file.path(conda_prefix, "lib")
  for (var in c("C_INCLUDE_PATH", "CPLUS_INCLUDE_PATH")) {
    cur <- Sys.getenv(var)
    if (!grepl(inc, cur, fixed = TRUE)) {
      new_val <- if (nchar(cur) > 0) paste(inc, cur, sep = ":") else inc
      do.call(Sys.setenv, setNames(list(new_val), var))
    }
  }
  cur_lib <- Sys.getenv("LIBRARY_PATH")
  if (!grepl(lib, cur_lib, fixed = TRUE))
    Sys.setenv(LIBRARY_PATH = if (nchar(cur_lib) > 0) paste(lib, cur_lib, sep = ":") else lib)
}

# bioconductor-motifmatchr / bioconductor-chromvar are not yet available on
# bioconda for R 4.5.x (Bioc 3.21), so install them from Bioconductor source.
if (!archr_done) {
  for (pkg in c("motifmatchr", "chromVAR")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("Installing", pkg, "via BiocManager...\n")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
  }

  remotes::install_github(
    "GreenleafLab/ArchR",
    ref   = "release_1.0.3",
    repos = BiocManager::repositories(),
    upgrade = "never"   # don't silently upgrade conda-managed packages
  )

  cat("\nVerifying installation...\n")
  library(ArchR)
  cat("ArchR", as.character(packageVersion("ArchR")), "installed successfully.\n")
}

# ── SeuratData + pbmcMultiModal reference ────────────────────────────────────
# SeuratData is not on conda or CRAN — install from GitHub.
# pbmcMultiModal is the Satija lab PBMC multimodal reference (Hao et al. 2021)
# used for celltype.l2 label transfer in the ALTRA pipeline.
cat("\n=== Installing SeuratData + pbmcMultiModal reference ===\n")

if (!requireNamespace("SeuratData", quietly = TRUE)) {
  cat("Installing SeuratData from GitHub...\n")
  remotes::install_github("satijalab/seurat-data", upgrade = "never", quiet = FALSE)
} else {
  cat("SeuratData", as.character(packageVersion("SeuratData")), "already installed.\n")
}

# Download pbmcref (Azimuth PBMC reference) if not already present.
# Hosted on seurat.nygenome.org — may be unreachable on HPC nodes without
# outbound internet. If the download fails, use the pre-converted RDS at
# data/ALTRA/pbmc_reference.rds (produced by bin/prep_reference.R).
if (!requireNamespace("pbmcref.SeuratData", quietly = TRUE)) {
  cat("Downloading pbmcref (Azimuth PBMC reference) via SeuratData...\n")
  library(SeuratData)
  tryCatch({
    SeuratData::InstallData("pbmcref")
    cat("pbmcref.SeuratData installed successfully.\n")
  }, error = function(e) {
    cat("NOTE: pbmcref download failed (likely no outbound network on this node).\n")
    cat("  Error:", conditionMessage(e), "\n")
    cat("  Alternative: use data/ALTRA/pbmc_reference.rds if it already exists,\n")
    cat("  or run bin/prep_reference.R in the 'seuratdisk' env to generate it.\n")
    cat("  SeuratData package is installed; only the reference data download was skipped.\n")
  })
} else {
  cat("pbmcref.SeuratData already installed.\n")
}

cat("\nLoad the PBMC reference with (if pbmcref.SeuratData download succeeded):\n")
cat("  library(SeuratData)\n")
cat("  reference <- SeuratData::LoadData('pbmcref')\n")
cat("\nOr load the pre-converted RDS directly:\n")
cat("  reference <- readRDS('data/ALTRA/pbmc_reference.rds')\n")
cat("\n=== postinstall_altra.R complete ===\n")
