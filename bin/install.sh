#!/usr/bin/env bash
# install.sh — one-command installer for the taiji-agent.
#
# Builds a conda env using composable PROFILES so users only install what
# their dataset actually needs. Then runs the R postinstall (SeuratDisk +
# MuDataSeurat from GitHub, only if the SC profile is present) and installs
# the Taiji binary. Idempotent: re-running on an existing install layers
# additional profiles on top without rebuilding from scratch.
#
# Profiles:
#   base   Python + xlsx + MACS3 + local taiji-agent package
#          Enables: detect-dataset-type, build-taiji-input, fetch-references,
#                   taiji-runner, workflow-log
#          ~500 MB,  ~5 min
#   sc     R + Seurat + Signac + Bioconductor (additive on top of base)
#          Enables: pseudobulk-construct
#          +3-4 GB, +15-30 min
#   dev    pytest + pytest-cov + ruff + mypy + ipython
#          Enables: development on the skills themselves
#          +500 MB, +3 min
#   full   = base + sc + dev (everything)
#
# Usage:
#   bash bin/install.sh                                    # base only (default)
#   bash bin/install.sh --system macos                     # base + macOS taiji binary
#   bash bin/install.sh --profile sc                       # base + sc
#   bash bin/install.sh --profile full                     # base + sc + dev
#   bash bin/install.sh --profile dev                      # base + dev
#   bash bin/install.sh --env-name my-env                  # non-default env name
#   bash bin/install.sh --skip-taiji                       # env-only
#   bash bin/install.sh --use-lockfile                     # install from conda-lock
#   bash bin/install.sh --solver mamba|micromamba|conda    # default: auto-detect
#                                                          # (prefers micromamba > mamba > conda)
#
# Exit codes:
#   0 ok
#   2 invalid args / missing solver
#   3 conda env create/update failed
#   4 postinstall.R failed
#   5 install-taiji.sh failed

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_BASE="${REPO_ROOT}/environment.base.yml"
ENV_SC="${REPO_ROOT}/environment.sc.yml"
ENV_DEV="${REPO_ROOT}/environment.dev.yml"
LOCKFILE="${REPO_ROOT}/conda-lock.linux-64.yml"

ENV_NAME="taiji-agent"
SYSTEM=""
PROFILE="base"
SKIP_TAIJI=0
SKIP_R=""    # tri-state: "" auto-decide, 1 force-skip, 0 force-run
USE_LOCKFILE=0
SOLVER=""    # empty means "auto-detect: prefer micromamba, then mamba, then conda"

usage() { sed -n '2,40p' "$0"; exit "${1:-0}"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --system)        SYSTEM="$2"; shift 2 ;;
    --profile)       PROFILE="$2"; shift 2 ;;
    --env-name)      ENV_NAME="$2"; shift 2 ;;
    --skip-taiji)    SKIP_TAIJI=1; shift ;;
    --skip-r)        SKIP_R=1; shift ;;
    --use-lockfile)  USE_LOCKFILE=1; shift ;;
    --solver)        SOLVER="$2"; shift 2 ;;
    -h|--help)       usage 0 ;;
    *) echo "unknown arg: $1" >&2; usage 2 ;;
  esac
done

# ---- Validate + expand profile ----
case "$PROFILE" in
  base)  PROFILES=(base) ;;
  sc)    PROFILES=(base sc) ;;
  dev)   PROFILES=(base dev) ;;
  full)  PROFILES=(base sc dev) ;;
  *) echo "unknown profile '$PROFILE' (expected: base|sc|dev|full)" >&2; exit 2 ;;
esac
echo "[install] profile=$PROFILE  (= ${PROFILES[*]})"

# Decide whether to run the R postinstall: only when the sc profile is in.
RUN_POSTINSTALL_R=0
for p in "${PROFILES[@]}"; do
  [[ "$p" == "sc" ]] && RUN_POSTINSTALL_R=1
done
# Honor explicit --skip-r override.
if [[ "$SKIP_R" == "1" ]]; then
  RUN_POSTINSTALL_R=0
fi

# ---- Resolve solver ----
# When the user didn't pass --solver explicitly, auto-detect in order:
# micromamba (fastest) -> mamba -> conda (slowest, but ~always available).
# When the user did pass --solver, honor it strictly (no fallback).
if [[ -z "$SOLVER" ]]; then
  for candidate in micromamba mamba conda; do
    if command -v "$candidate" >/dev/null 2>&1; then
      SOLVER="$candidate"
      break
    fi
  done
  if [[ -z "$SOLVER" ]]; then
    echo "[install] no conda solver found on PATH (tried: micromamba, mamba, conda)" >&2
    echo "[install] install one of:" >&2
    echo "    micromamba (fastest):  https://mamba.readthedocs.io/en/latest/micromamba-installation.html" >&2
    echo "    miniforge/mamba:       https://github.com/conda-forge/miniforge" >&2
    echo "    miniconda:             https://docs.conda.io/en/latest/miniconda.html" >&2
    exit 2
  fi
  echo "[install] auto-detected solver: $SOLVER ($(command -v "$SOLVER"))"
else
  if ! command -v "$SOLVER" >/dev/null 2>&1; then
    echo "[install] solver '$SOLVER' (passed via --solver) not on PATH." >&2
    echo "[install] either install it, or omit --solver to auto-detect from " >&2
    echo "[install] micromamba | mamba | conda." >&2
    exit 2
  fi
  echo "[install] using solver: $SOLVER ($(command -v "$SOLVER"))"
fi

# ---- Map profile -> env-file path ----
profile_file() {
  case "$1" in
    base) echo "$ENV_BASE" ;;
    sc)   echo "$ENV_SC" ;;
    dev)  echo "$ENV_DEV" ;;
    *)    echo "" ;;
  esac
}

# Verify profile files exist before doing anything.
for p in "${PROFILES[@]}"; do
  f=$(profile_file "$p")
  if [[ ! -f "$f" ]]; then
    echo "[install] profile file missing: $f" >&2
    exit 2
  fi
done

# ---- Lockfile path takes priority if --use-lockfile is set ----
if [[ "$USE_LOCKFILE" -eq 1 ]]; then
  if [[ ! -f "$LOCKFILE" ]]; then
    echo "[install] --use-lockfile requested but $LOCKFILE missing." >&2
    echo "[install] generate it with: conda-lock --file environment.yml --platform linux-64" >&2
    exit 2
  fi
  echo "[install] env source: $LOCKFILE (lockfile, full env)"
  if "$SOLVER" env list 2>/dev/null | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    "$SOLVER" remove -y -n "$ENV_NAME" --all
  fi
  "$SOLVER" create -y -n "$ENV_NAME" -f "$LOCKFILE" \
    || { echo "[install] env create from lockfile failed" >&2; exit 3; }
else
  # ---- Profile-by-profile install (additive) ----
  ENV_EXISTS=0
  if "$SOLVER" env list 2>/dev/null | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    ENV_EXISTS=1
  fi

  for p in "${PROFILES[@]}"; do
    f=$(profile_file "$p")
    if [[ "$ENV_EXISTS" -eq 0 ]]; then
      echo "[install] creating env '$ENV_NAME' from $f"
      "$SOLVER" create -y -n "$ENV_NAME" -f "$f" \
        || { echo "[install] env create failed (profile $p)" >&2; exit 3; }
      ENV_EXISTS=1
    else
      echo "[install] layering profile '$p' onto existing env '$ENV_NAME' ($f)"
      "$SOLVER" install -y -n "$ENV_NAME" -f "$f" \
        || { echo "[install] env install failed (profile $p)" >&2; exit 3; }
    fi
  done
fi

# ---- Activate-then-run helper ----
# `conda run` historically buffers stdout/stderr; --no-capture-output fixes
# that. micromamba/mamba run don't have the issue.
if [[ "$SOLVER" == "conda" ]]; then
  RUN=("$SOLVER" run --no-capture-output -n "$ENV_NAME")
else
  RUN=("$SOLVER" run -n "$ENV_NAME")
fi

# ---- R postinstall (only when sc profile is in scope) ----
if [[ "$RUN_POSTINSTALL_R" -eq 1 ]]; then
  echo "[install] running R postinstall (SeuratDisk + MuDataSeurat) ..."
  if ! "${RUN[@]}" Rscript "${REPO_ROOT}/bin/postinstall.R"; then
    echo "[install] postinstall.R failed; investigate before proceeding." >&2
    exit 4
  fi
else
  if [[ "$SKIP_R" == "1" ]]; then
    echo "[install] --skip-r set; skipping GitHub R-package install."
  else
    echo "[install] sc profile not in scope; skipping R postinstall."
  fi
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

# ---- Final summary ----
TAIJI_PATH="(skipped)"
[[ "$SKIP_TAIJI" -eq 0 ]] && TAIJI_PATH="${REPO_ROOT}/binaries/taiji"

cat <<EOF

[install] SUCCESS
  env name:       $ENV_NAME
  profile:        $PROFILE  (=${PROFILES[*]})
  Taiji binary:   $TAIJI_PATH

Activate with:
  $SOLVER activate $ENV_NAME

Verify everything is wired up:
  bash bin/doctor.sh --profile $PROFILE
EOF

# ---- Profile-upgrade hint ----
if [[ "$PROFILE" == "base" ]]; then
  cat <<EOF

NOTE: you installed the 'base' profile (bulk-only).
If you later need single-cell skills (pseudobulk-construct), upgrade with:
  bash bin/install.sh --profile sc
which adds R + Seurat/Signac on top of the existing env (no full reinstall).
EOF
fi
