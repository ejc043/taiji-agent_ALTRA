#!/usr/bin/env bash
# install.sh — one-command installer for the taiji-agent.
#
# Builds the conda env, runs the R postinstall (SeuratDisk + MuDataSeurat from
# GitHub), and installs the Taiji binary for the chosen system. Idempotent:
# re-running on an existing install updates whatever needs updating and is a
# no-op for what's already in place.
#
# Usage:
#   bash bin/install.sh                              # full install, auto-detect OS
#   bash bin/install.sh --system centos              # explicit Taiji binary target
#   bash bin/install.sh --env-name my-taiji-env      # non-default env name
#   bash bin/install.sh --skip-taiji                 # env-only, no Taiji binary
#   bash bin/install.sh --skip-r                     # env + Taiji, skip postinstall.R
#   bash bin/install.sh --use-lockfile               # install from conda-lock.linux-64.yml
#   bash bin/install.sh --solver mamba|micromamba|conda    # default: micromamba
#
# Exit codes:
#   0 ok
#   2 invalid args / missing solver binary
#   3 conda env create/update failed
#   4 postinstall.R failed
#   5 install-taiji.sh failed

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_FILE="${REPO_ROOT}/environment.yml"
LOCKFILE="${REPO_ROOT}/conda-lock.linux-64.yml"

ENV_NAME="taiji-agent"
SYSTEM=""
SKIP_TAIJI=0
SKIP_R=0
USE_LOCKFILE=0
SOLVER="micromamba"

usage() { sed -n '2,17p' "$0"; exit "${1:-0}"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --system)        SYSTEM="$2"; shift 2 ;;
    --env-name)      ENV_NAME="$2"; shift 2 ;;
    --skip-taiji)    SKIP_TAIJI=1; shift ;;
    --skip-r)        SKIP_R=1; shift ;;
    --use-lockfile)  USE_LOCKFILE=1; shift ;;
    --solver)        SOLVER="$2"; shift 2 ;;
    -h|--help)       usage 0 ;;
    *) echo "unknown arg: $1" >&2; usage 2 ;;
  esac
done

# ---- Resolve solver binary ----
if ! command -v "$SOLVER" >/dev/null 2>&1; then
  echo "[install] solver '$SOLVER' not on PATH." >&2
  echo "[install] install micromamba: https://mamba.readthedocs.io/en/latest/micromamba-installation.html" >&2
  echo "[install] or pass --solver mamba|conda if you have one of those." >&2
  exit 2
fi
echo "[install] using solver: $SOLVER ($(command -v "$SOLVER"))"

# ---- Pick the source of truth: lockfile or environment.yml ----
SRC_FILE="$ENV_FILE"
if [[ "$USE_LOCKFILE" -eq 1 ]]; then
  if [[ ! -f "$LOCKFILE" ]]; then
    echo "[install] --use-lockfile requested but $LOCKFILE missing." >&2
    echo "[install] generate it with: conda-lock --file environment.yml --platform linux-64" >&2
    exit 2
  fi
  SRC_FILE="$LOCKFILE"
fi
echo "[install] env source: $SRC_FILE"

# ---- Create or update the env ----
# `micromamba env list` is the portable way to check existence across solvers.
if "$SOLVER" env list 2>/dev/null | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "[install] env '$ENV_NAME' exists — updating in place."
  # `update` for environment.yml; for lockfile-based installs, recreating is
  # the canonical path because lockfiles are immutable snapshots.
  if [[ "$USE_LOCKFILE" -eq 1 ]]; then
    "$SOLVER" remove -y -n "$ENV_NAME" --all
    "$SOLVER" create  -y -n "$ENV_NAME" -f "$SRC_FILE"
  else
    "$SOLVER" install -y -n "$ENV_NAME" -f "$SRC_FILE"
  fi
else
  echo "[install] creating env '$ENV_NAME' from $SRC_FILE"
  "$SOLVER" create -y -n "$ENV_NAME" -f "$SRC_FILE"
fi || { echo "[install] env create/update failed" >&2; exit 3; }

# ---- Activate-then-run helper. micromamba's `run` subcommand is the most
# portable way to execute a command inside the env without sourcing
# activate scripts (which differ across shells).
RUN=("$SOLVER" run -n "$ENV_NAME")

# ---- R postinstall: SeuratDisk + MuDataSeurat from GitHub ----
if [[ "$SKIP_R" -eq 0 ]]; then
  echo "[install] running R postinstall (SeuratDisk + MuDataSeurat) ..."
  if ! "${RUN[@]}" Rscript "${REPO_ROOT}/bin/postinstall.R"; then
    echo "[install] postinstall.R failed; investigate before proceeding." >&2
    exit 4
  fi
else
  echo "[install] --skip-r set; skipping GitHub R-package install."
fi

# ---- Taiji binary ----
if [[ "$SKIP_TAIJI" -eq 0 ]]; then
  TAIJI_ARGS=()
  [[ -n "$SYSTEM" ]] && TAIJI_ARGS+=(--system "$SYSTEM")
  if ! bash "${REPO_ROOT}/bin/install-taiji.sh" "${TAIJI_ARGS[@]}"; then
    echo "[install] install-taiji.sh failed" >&2
    exit 5
  fi
else
  echo "[install] --skip-taiji set; not downloading the Taiji binary."
fi

cat <<EOF

[install] SUCCESS
  env name:       $ENV_NAME
  source:         $SRC_FILE
  Taiji binary:   $( [[ "$SKIP_TAIJI" -eq 0 ]] && echo "${REPO_ROOT}/binaries/taiji" || echo "(skipped)" )

Activate with:
  $SOLVER activate $ENV_NAME

Verify everything is wired up:
  bash bin/doctor.sh
EOF
