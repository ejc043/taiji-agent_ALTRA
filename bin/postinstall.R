#!/usr/bin/env Rscript
# postinstall.R — install GitHub-only R packages after the conda env is built.
#
# Why this isn't in environment.yml: SeuratDisk and MuDataSeurat are not
# distributed via bioconda or conda-forge. The official install path is
# remotes::install_github, which is straightforward but can't be expressed
# in a conda recipe.
#
# Idempotency: this script checks whether each package already loads cleanly
# before reinstalling. Re-runs are a no-op when packages are present.
#
# Usage:
#   Rscript bin/postinstall.R
#   Rscript bin/postinstall.R --force        # reinstall even if already loadable

args <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args

suppressPackageStartupMessages({
  if (!requireNamespace("remotes", quietly = TRUE)) {
    stop("remotes package missing — was the conda env built? ",
         "Run `micromamba install -n taiji-agent r-remotes` first.")
  }
})

# Each entry: (R package name, GitHub repo). Add new GitHub-only deps here.
github_packages <- list(
  list(pkg = "SeuratDisk",   repo = "mojaveazure/seurat-disk"),
  list(pkg = "MuDataSeurat", repo = "PMBio/MuDataSeurat")
)

install_one <- function(pkg, repo) {
  loadable <- requireNamespace(pkg, quietly = TRUE)
  if (loadable && !force) {
    message(sprintf("[postinstall] %s already installed (%s) — skipping.",
                    pkg, packageVersion(pkg)))
    return(invisible(TRUE))
  }
  message(sprintf("[postinstall] installing %s from github:%s ...", pkg, repo))
  tryCatch(
    {
      remotes::install_github(repo, upgrade = "never", quiet = FALSE,
                              dependencies = TRUE)
      # Verify it loads after install.
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("%s installed but does not load. Inspect the build log above.", pkg))
      }
      message(sprintf("[postinstall] OK: %s %s", pkg, packageVersion(pkg)))
      invisible(TRUE)
    },
    error = function(e) {
      message(sprintf("[postinstall] FAILED: %s — %s", pkg, conditionMessage(e)))
      invisible(FALSE)
    }
  )
}

results <- vapply(github_packages,
                  function(x) install_one(x$pkg, x$repo),
                  logical(1))

if (any(!results)) {
  cat("\n[postinstall] one or more packages failed. See messages above.\n",
      file = stderr())
  quit(status = 1)
}

cat("\n[postinstall] all GitHub R packages installed successfully.\n")
