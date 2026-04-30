#!/usr/bin/env bash
# install.sh — one-command installer for the taiji-agent.
#
# Builds a conda env using composable PROFILES so users only install what
# their dataset actually needs. Then runs the R postinstall (MuDataSeurat
# from GitHub, only if the SC profile is present) and installs the Taiji
# binary. Idempotent: re-running on an existing install layers additional
# profiles on top without rebuilding from scratch.
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
#
# Usage:
#   bash bin/install.sh                                    # base only (default)
#   bash bin/install.sh --system macos                     # base + macOS taiji binary
#   bash bin/install.sh --profile sc                       # base + sc
#   bash bin/install.sh --profile dev                      # base + dev
#   bash bin/install.sh --env-name my-env                  # non-default env name
#   bash bin/install.sh --env-prefix /abs/path/envs/foo    # explicit env path
#                                                          # (overrides auto-detect;
#                                                          #  use when the env is at
#                                                          #  a non-standard location
#                                                          #  e.g. /stg3/.../envs/...)
#   bash bin/install.sh --skip-taiji                       # env-only
#   bash bin/install.sh --use-lockfile                     # install from conda-lock
#   bash bin/install.sh --solver mamba|micromamba|conda    # default: auto-detect
#                                                          # (prefers micromamba > mamba > conda)
#   bash bin/install.sh --system macos --fetch-references --genome hg38
#                                                          # install + stage reference data
#                                                          # in one step (no separate activate)
#
# Exit codes:
#   0 ok
#   2 invalid args / missing solver
#   3 conda env create/update failed
#   4 postinstall.R failed
#   5 install-taiji.sh failed
#   6 fetch-references failed

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_BASE="${REPO_ROOT}/environment.base.yml"
ENV_SC="${REPO_ROOT}/environment.sc.yml"
ENV_DEV="${REPO_ROOT}/environment.dev.yml"
LOCKFILE="${REPO_ROOT}/conda-lock.linux-64.yml"

ENV_NAME="taiji-agent"
ENV_PREFIX=""    # explicit absolute path override (--env-prefix); empty = auto-detect
SYSTEM=""
PROFILE="base"
SKIP_TAIJI=0
SKIP_R=""    # tri-state: "" auto-decide, 1 force-skip, 0 force-run
USE_LOCKFILE=0
SOLVER=""    # empty means "auto-detect: prefer micromamba, then mamba, then conda"
FETCH_REFS=0
FETCH_GENOME=""
FETCH_OUTPUT=""    # empty = default (dependencies_data/ inside repo root)

usage() { sed -n '2,40p' "$0"; exit "${1:-0}"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --system)        SYSTEM="$2"; shift 2 ;;
    --profile)       PROFILE="$2"; shift 2 ;;
    --env-name)      ENV_NAME="$2"; shift 2 ;;
    --env-prefix)    ENV_PREFIX="$2"; shift 2 ;;
    --skip-taiji)         SKIP_TAIJI=1; shift ;;
    --skip-r)             SKIP_R=1; shift ;;
    --use-lockfile)       USE_LOCKFILE=1; shift ;;
    --solver)             SOLVER="$2"; shift 2 ;;
    --fetch-references)   FETCH_REFS=1; shift ;;
    --genome)             FETCH_GENOME="$2"; shift 2 ;;
    --fetch-output)       FETCH_OUTPUT="$2"; shift 2 ;;
    -h|--help)            usage 0 ;;
    *) echo "unknown arg: $1" >&2; usage 2 ;;
  esac
done

# ---- Validate + expand profile ----
case "$PROFILE" in
  base)  PROFILES=(base) ;;
  sc)    PROFILES=(base sc) ;;
  dev)   PROFILES=(base dev) ;;
  *) echo "unknown profile '$PROFILE' (expected: base|sc|dev)" >&2; exit 2 ;;
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

# ---- Profile-cache helpers --------------------------------------------------
#
# We track which profiles are already installed in the conda env via a
# marker file inside the env's directory. Each line is:
#     <profile_name>  <12-char sha256 hash of the env-file when installed>
#
# Re-running `bin/install.sh --profile X` is then a no-op when the marker
# already records X with the current hash. Editing environment.<X>.yml
# changes the hash and forces a re-install.

# Compute the conda env's filesystem path; empty if env doesn't exist.
#
# Tries multiple discovery mechanisms in order, because micromamba/conda
# only see envs under their CURRENT root prefix — but on shared
# infrastructure (e.g. UCSD SLURM with envs under /stg3/) the active
# shell's MAMBA_ROOT_PREFIX may differ from where the env actually lives:
#
#   1. --env-prefix flag (user-specified absolute path)
#   2. $CONDA_PREFIX matches ENV_NAME (env is currently activated)
#   3. <solver> env list (works when MAMBA_ROOT_PREFIX matches the install)
#   4. Probe common candidate paths (HOME-relative, $MAMBA_ROOT_PREFIX,
#      site-specific /stg3/, /opt/, etc.)
get_env_path() {
  # 1. Explicit override
  if [[ -n "${ENV_PREFIX:-}" ]]; then
    [[ -d "$ENV_PREFIX" ]] && { echo "$ENV_PREFIX"; return 0; }
    echo "[install] WARN: --env-prefix '$ENV_PREFIX' does not exist" >&2
    return 1
  fi

  # 2. Currently activated env
  if [[ -n "${CONDA_PREFIX:-}" && "$(basename "$CONDA_PREFIX" 2>/dev/null)" == "$ENV_NAME" ]]; then
    echo "$CONDA_PREFIX"
    return 0
  fi

  # 3. Solver's own listing (only sees its current root prefix)
  local from_solver
  from_solver=$("$SOLVER" env list 2>/dev/null \
    | awk -v e="$ENV_NAME" '$1==e {print $NF}' \
    | head -1)
  if [[ -n "$from_solver" && -d "$from_solver" ]]; then
    echo "$from_solver"
    return 0
  fi

  # 4. Probe common candidate locations.
  local user="${USER:-$(id -un 2>/dev/null)}"
  local candidates=(
    "${MAMBA_ROOT_PREFIX:-}/envs/$ENV_NAME"
    "${CONDA_ENVS_PATH:-}/$ENV_NAME"
    "$HOME/.local/share/mamba/envs/$ENV_NAME"
    "$HOME/micromamba/envs/$ENV_NAME"
    "$HOME/miniforge3/envs/$ENV_NAME"
    "$HOME/miniconda3/envs/$ENV_NAME"
    "$HOME/anaconda3/envs/$ENV_NAME"
    "$HOME/opt/miniconda3/envs/$ENV_NAME"
    "$HOME/opt/anaconda3/envs/$ENV_NAME"
    "/opt/conda/envs/$ENV_NAME"
    "/opt/miniconda3/envs/$ENV_NAME"
    "/opt/anaconda3/envs/$ENV_NAME"
    # UCSD-specific shared scratch (Eunice's setup)
    "/stg3/data1/$user/.local/share/mamba/envs/$ENV_NAME"
    "/stg3/data1/$user/micromamba/envs/$ENV_NAME"
    "/stg3/data1/$user/miniconda3/envs/$ENV_NAME"
  )
  local c
  for c in "${candidates[@]}"; do
    [[ -n "$c" && "$c" != "/envs/$ENV_NAME" && -d "$c" ]] && { echo "$c"; return 0; }
  done

  # Not found anywhere
  return 1
}

# Print where we looked, for debugging "env not detected" surprises.
debug_env_search() {
  echo "[install] env-detection debug — checking these locations for '$ENV_NAME':"
  echo "  --env-prefix:                ${ENV_PREFIX:-(unset)}"
  echo "  \$CONDA_PREFIX (activated):    ${CONDA_PREFIX:-(unset)}"
  echo "  $SOLVER env list:"
  "$SOLVER" env list 2>/dev/null | sed 's/^/    /'
  echo "  candidate paths (probed in order):"
  local user="${USER:-$(id -un 2>/dev/null)}"
  for c in \
    "${MAMBA_ROOT_PREFIX:-}/envs/$ENV_NAME" \
    "${CONDA_ENVS_PATH:-}/$ENV_NAME" \
    "$HOME/.local/share/mamba/envs/$ENV_NAME" \
    "$HOME/micromamba/envs/$ENV_NAME" \
    "$HOME/miniforge3/envs/$ENV_NAME" \
    "$HOME/miniconda3/envs/$ENV_NAME" \
    "$HOME/anaconda3/envs/$ENV_NAME" \
    "/opt/conda/envs/$ENV_NAME" \
    "/stg3/data1/$user/.local/share/mamba/envs/$ENV_NAME" \
    "/stg3/data1/$user/micromamba/envs/$ENV_NAME"; do
    if [[ -n "$c" && "$c" != "/envs/$ENV_NAME" ]]; then
      [[ -d "$c" ]] && echo "    [FOUND] $c" || echo "    [miss]  $c"
    fi
  done
}

# Short SHA-256 of a file (12 hex chars). Uses python3 for portability
# (sha256sum/shasum availability differs across macOS vs Linux).
file_hash() {
  python3 - "$1" <<'PY' 2>/dev/null || echo "nohash"
import hashlib, sys
print(hashlib.sha256(open(sys.argv[1], "rb").read()).hexdigest()[:12])
PY
}

# Read the marker file (if present) into a global associative-array-like
# string-map (bash 3 friendly): "$PROFILE_STATE" lines = "<profile> <hash>".
read_profile_state() {
  local env_path="$1"
  local marker="$env_path/.taiji-agent-profiles"
  [[ -f "$marker" ]] && cat "$marker" || echo ""
}

# Append/update a profile entry in the marker.
write_profile_state() {
  local env_path="$1" profile="$2" hash="$3"
  local marker="$env_path/.taiji-agent-profiles"
  # Drop any existing line for this profile, then append the new one.
  if [[ -f "$marker" ]]; then
    grep -v "^${profile} " "$marker" > "${marker}.tmp" 2>/dev/null || true
    mv "${marker}.tmp" "$marker"
  fi
  echo "${profile} ${hash}" >> "$marker"
}

# True if the marker says <profile> is installed AND its hash matches the
# current env file's hash.
profile_already_installed() {
  local env_path="$1" profile="$2"
  local f cur_hash recorded_hash state
  [[ -z "$env_path" || ! -d "$env_path" ]] && return 1
  f=$(profile_file "$profile")
  [[ ! -f "$f" ]] && return 1
  cur_hash=$(file_hash "$f")
  state=$(read_profile_state "$env_path")
  recorded_hash=$(echo "$state" | awk -v p="$profile" '$1==p {print $2}' | head -1)
  [[ -n "$recorded_hash" && "$recorded_hash" == "$cur_hash" ]]
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
  EXISTING_PATH=$(get_env_path 2>/dev/null || true)
  if [[ -n "$EXISTING_PATH" && -d "$EXISTING_PATH" ]]; then
    "$SOLVER" remove -y --prefix "$EXISTING_PATH" --all
  fi
  "$SOLVER" create -y -n "$ENV_NAME" -f "$LOCKFILE" \
    || { echo "[install] env create from lockfile failed" >&2; exit 3; }
else
  # ---- Profile-by-profile install (additive + cache-aware) ----
  ENV_EXISTS=0
  ENV_PATH=""
  # get_env_path probes multiple discovery paths (--env-prefix, $CONDA_PREFIX,
  # `<solver> env list`, and common candidate dirs incl. UCSD's /stg3/...)
  # so we find the env even when MAMBA_ROOT_PREFIX in the current shell
  # doesn't point at where the env was originally created.
  if ENV_PATH=$(get_env_path); then
    ENV_EXISTS=1
    echo "[install] env '$ENV_NAME' already exists at: $ENV_PATH"
  else
    echo "[install] env '$ENV_NAME' not found in any standard location"
    echo "[install]   (run with TAIJI_DEBUG_ENV=1 or pass --env-prefix to override)"
  fi
  if [[ -n "${TAIJI_DEBUG_ENV:-}" ]]; then
    debug_env_search
  fi

  for p in "${PROFILES[@]}"; do
    f=$(profile_file "$p")

    # ---- Cache-hit short-circuit: if the env exists AND the marker says
    # this profile was installed AND the env file's hash hasn't changed
    # since then, skip the conda call entirely.
    if [[ "$ENV_EXISTS" -eq 1 ]] && profile_already_installed "$ENV_PATH" "$p"; then
      echo "[install] profile '$p' already installed (env-file hash matches); skipping."
      continue
    fi

    if [[ "$ENV_EXISTS" -eq 0 ]]; then
      echo "[install] creating env '$ENV_NAME' from $f"
      "$SOLVER" create -y -n "$ENV_NAME" -f "$f" \
        || { echo "[install] env create failed (profile $p)" >&2; exit 3; }
      ENV_EXISTS=1
      ENV_PATH=$(get_env_path)
    else
      # Use --prefix when we have an absolute path; the solver's -n lookup
      # only sees envs under its current root prefix, but get_env_path
      # may have found the env at a different root (e.g. /stg3/...).
      echo "[install] layering profile '$p' onto existing env '$ENV_PATH' ($f)"
      if [[ -n "$ENV_PATH" && -d "$ENV_PATH" ]]; then
        "$SOLVER" install -y --prefix "$ENV_PATH" -f "$f" \
          || { echo "[install] env install failed (profile $p)" >&2; exit 3; }
      else
        "$SOLVER" install -y -n "$ENV_NAME" -f "$f" \
          || { echo "[install] env install failed (profile $p)" >&2; exit 3; }
      fi
    fi

    # Record this profile in the marker with the current env-file hash so
    # subsequent runs short-circuit.
    if [[ -n "$ENV_PATH" && -d "$ENV_PATH" ]]; then
      write_profile_state "$ENV_PATH" "$p" "$(file_hash "$f")"
    fi
  done
fi

# ---- Activate-then-run helper ----
# `conda run` historically buffers stdout/stderr; --no-capture-output fixes
# that. micromamba/mamba run don't have the issue.
# Prefer --prefix <path> over -n <name> when we have an absolute env path —
# avoids "env not found" on shells where MAMBA_ROOT_PREFIX doesn't point
# at the env's actual location.
ENV_SELECTOR=()
if [[ -n "$ENV_PATH" && -d "$ENV_PATH" ]]; then
  ENV_SELECTOR=(--prefix "$ENV_PATH")
else
  ENV_SELECTOR=(-n "$ENV_NAME")
fi
if [[ "$SOLVER" == "conda" ]]; then
  RUN=("$SOLVER" run --no-capture-output "${ENV_SELECTOR[@]}")
else
  RUN=("$SOLVER" run "${ENV_SELECTOR[@]}")
fi

# ---- R postinstall (only when sc profile is in scope) ----
if [[ "$RUN_POSTINSTALL_R" -eq 1 ]]; then
  echo "[install] running R postinstall (MuDataSeurat) ..."
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

# ---- Stage reference data (only when --fetch-references is set) ----
if [[ "$FETCH_REFS" -eq 1 ]]; then
  if [[ -z "$FETCH_GENOME" ]]; then
    echo "[install] --fetch-references requires --genome <build>" >&2
    echo "[install]   supported builds: hg38, hg19, mm10" >&2
    echo "[install]   e.g. bash bin/install.sh --fetch-references --genome hg38" >&2
    exit 2
  fi
  FETCH_OUT="${FETCH_OUTPUT:-${REPO_ROOT}/dependencies_data}"
  FETCH_PY="${REPO_ROOT}/skills/fetch-references/scripts/fetch.py"
  echo "[install] staging reference data for $FETCH_GENOME into $FETCH_OUT ..."
  if ! "${RUN[@]}" python "$FETCH_PY" \
      --genome "$FETCH_GENOME" \
      --output "$FETCH_OUT" \
      --update-genomes-yml; then
    echo "[install] fetch-references failed" >&2
    exit 6
  fi
else
  echo "[install] skipping reference data (pass --fetch-references --genome <build> to stage in one step)."
fi

# ---- Register env's parent dir in solver's envs_dirs ----
# So `<solver> activate <name>` works in any future shell, regardless of
# whatever MAMBA_ROOT_PREFIX is set to. By default, micromamba/mamba/conda
# only resolve env names under their current root prefix's envs/ subdir;
# on clusters where the env lives under a non-default prefix (here:
# /stg3/...) the user has to type the full path every time.
# Appending to envs_dirs writes a persistent entry into the solver's
# config (~/.condarc / ~/.mambarc) so name lookup always finds this env.
register_envs_dir() {
  local env_path="$1"
  [[ -z "$env_path" || ! -d "$env_path" ]] && return 0
  local parent
  parent=$(dirname "$env_path")
  # Already in envs_dirs? skip to avoid duplicate entries on re-runs.
  if "$SOLVER" config list 2>/dev/null | grep -Fq -- "- $parent"; then
    echo "[install] envs_dirs already registers '$parent'; skipping."
    return 0
  fi
  # Command syntax differs across solvers:
  #   micromamba: `config append envs_dirs <path>`
  #   mamba/conda: `config --append envs_dirs <path>`
  local rc
  case "$SOLVER" in
    micromamba)
      "$SOLVER" config append envs_dirs "$parent" >/dev/null 2>&1; rc=$?
      ;;
    mamba|conda)
      "$SOLVER" config --append envs_dirs "$parent" >/dev/null 2>&1; rc=$?
      ;;
    *) rc=99 ;;
  esac
  if [[ "$rc" -eq 0 ]]; then
    echo "[install] registered '$parent' in $SOLVER envs_dirs"
    echo "[install]   ($SOLVER activate $ENV_NAME will now work in any shell)"
  else
    echo "[install] NOTE: could not auto-register envs_dirs via $SOLVER config." >&2
    echo "[install]   Activation by name may fail in fresh shells; use the" >&2
    echo "[install]   full-prefix form instead: $SOLVER activate $env_path" >&2
  fi
}
register_envs_dir "${ENV_PATH:-}"

# ---- Final summary ----
TAIJI_PATH="(skipped)"
[[ "$SKIP_TAIJI" -eq 0 ]] && TAIJI_PATH="${REPO_ROOT}/binaries/taiji"

REF_DATA_LINE="(not staged — rerun with --fetch-references --genome <build> to stage)"
if [[ "$FETCH_REFS" -eq 1 ]]; then
  REF_DATA_LINE="${FETCH_OUTPUT:-${REPO_ROOT}/dependencies_data}/$FETCH_GENOME/"
fi

cat <<EOF

[install] SUCCESS
  env name:       $ENV_NAME
  env path:       ${ENV_PATH:-(unknown)}
  profile:        $PROFILE  (=${PROFILES[*]})
  Taiji binary:   $TAIJI_PATH
  reference data: $REF_DATA_LINE

Activate with (in this order of preference):
  $SOLVER activate $ENV_NAME            # by-name (works in any shell after the envs_dirs registration above)
  $SOLVER activate ${ENV_PATH:-<env_path>}   # by-prefix (always works, even before envs_dirs takes effect)

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
