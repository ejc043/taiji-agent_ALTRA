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

if (requireNamespace("ArchR", quietly = TRUE)) {
  cat("ArchR is already installed:", as.character(packageVersion("ArchR")), "\n")
  cat("Nothing to do.\n")
  quit(status = 0)
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  stop("remotes package not found. Is taiji-agent-altra activated?\n",
       "Run: micromamba activate taiji-agent-altra")
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
