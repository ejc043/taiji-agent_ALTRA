#!/usr/bin/env bash
# run-taiji.sh — top-level Taiji executor.
#
# Thin shell wrapper around skills/taiji-runner. Picks the right Taiji
# binary for this OS+arch, regenerates the xlsx with absolute paths to
# THIS machine, runs preflight checks, starts a workflow-log run, then
# delegates per-sample TSV/config materialization + per-sample Taiji
# invocation to skills/taiji-runner/scripts/run_taiji.py.
#
# Usage:
#   bash bin/run-taiji.sh                                   # default run_dir
#   bash bin/run-taiji.sh runs/<name>                       # other run dir
#   bash bin/run-taiji.sh runs/<name> --threads 8           # parallelism
#   bash bin/run-taiji.sh runs/<name> --parallel 2          # 2 samples concurrent
#   bash bin/run-taiji.sh runs/<name> --samples RA_11       # subset
#   bash bin/run-taiji.sh runs/<name> --prepare-only        # no taiji invoked
#   bash bin/run-taiji.sh runs/<name> --skip-preflight      # for re-runs
#   bash bin/run-taiji.sh runs/<name> --dry-run             # plan only
#   bash bin/run-taiji.sh runs/<name> --binary <abs path>   # override binary
#   bash bin/run-taiji.sh runs/<name> --strict-sandbox      # require OS-level sandbox
#                                                           # (bwrap on Linux,
#                                                           #  sandbox-exec on macOS)
#   bash bin/run-taiji.sh runs/<name> --no-sandbox          # bypass the wrapper entirely
#                                                           # (debugging only)
#
# Executor flags (forwarded to taiji-runner):
#   --executor {auto|local|slurm}    auto = slurm if sbatch on PATH else local
#   --max-parallel N                 SLURM array concurrency cap (default 10)
#   --slurm-mem-gb N                 SLURM per-task memory (default 30)
#   --slurm-time HH:MM:SS            SLURM per-task time limit (default 24:00:00)
#   --no-wait                        SLURM only: submit + return; don't poll squeue
#
# Env vars: THREADS=N (default 4), PARALLEL=N (default 1 local / 10 slurm),
#           TAIJI_STRICT_SANDBOX=1 (same as --strict-sandbox)
#
# Exit codes: 0 ok, 2 bad args, 3 missing dep, 4 preflight failed,
#             5 one or more samples failed.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# ---- Args ----------------------------------------------------------------

RUN_DIR="${1:-$REPO_ROOT/runs/RA_OA_chr22_demo}"
shift || true

SKIP_PREFLIGHT=0
DRY_RUN=0
PREPARE_ONLY=0
BINARY_OVERRIDE=""
SAMPLES=""
EXTRA_RUNNER_ARGS=()
STRICT_SANDBOX=0
NO_SANDBOX=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --skip-preflight)     SKIP_PREFLIGHT=1; shift ;;
    --dry-run)            DRY_RUN=1; shift ;;
    --prepare-only)       PREPARE_ONLY=1; shift ;;
    --binary)             BINARY_OVERRIDE="$2"; shift 2 ;;
    --samples)            SAMPLES="$2"; shift 2 ;;
    --threads)            EXTRA_RUNNER_ARGS+=(--threads "$2"); shift 2 ;;
    --parallel)           EXTRA_RUNNER_ARGS+=(--parallel "$2"); shift 2 ;;
    --executor)           EXTRA_RUNNER_ARGS+=(--executor "$2"); shift 2 ;;
    --max-parallel)       EXTRA_RUNNER_ARGS+=(--max-parallel "$2"); shift 2 ;;
    --slurm-mem-gb)       EXTRA_RUNNER_ARGS+=(--slurm-mem-gb "$2"); shift 2 ;;
    --slurm-time)         EXTRA_RUNNER_ARGS+=(--slurm-time "$2"); shift 2 ;;
    --no-wait)            EXTRA_RUNNER_ARGS+=(--no-wait); shift ;;
    --continue-on-error)  EXTRA_RUNNER_ARGS+=(--continue-on-error); shift ;;
    --strict-sandbox)     STRICT_SANDBOX=1; shift ;;
    --no-sandbox)         NO_SANDBOX=1; shift ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

THREADS="${THREADS:-4}"
PARALLEL="${PARALLEL:-1}"
# If user didn't pass --threads/--parallel explicitly, use env vars.
[[ " ${EXTRA_RUNNER_ARGS[*]:-} " == *" --threads "* ]] || \
    EXTRA_RUNNER_ARGS+=(--threads "$THREADS")
[[ " ${EXTRA_RUNNER_ARGS[*]:-} " == *" --parallel "* ]] || \
    EXTRA_RUNNER_ARGS+=(--parallel "$PARALLEL")

RUN_DIR="$(cd "$RUN_DIR" && pwd)"
[[ -d "$RUN_DIR" ]] || { echo "ERROR: run dir does not exist: $RUN_DIR" >&2; exit 2; }

# Resolve samples.csv + taiji_config.template.yml. Canonical location is
# <run_dir>/code/; flat layout (<run_dir>/...) is accepted for back-compat.
if   [[ -f "$RUN_DIR/code/samples.csv" ]]; then SAMPLES_CSV="$RUN_DIR/code/samples.csv"
elif [[ -f "$RUN_DIR/samples.csv"      ]]; then SAMPLES_CSV="$RUN_DIR/samples.csv"
else echo "ERROR: samples.csv missing (looked in $RUN_DIR/code/ and $RUN_DIR/)" >&2; exit 2
fi
if   [[ -f "$RUN_DIR/code/taiji_config.template.yml" ]]; then TEMPLATE_YML="$RUN_DIR/code/taiji_config.template.yml"
elif [[ -f "$RUN_DIR/taiji_config.template.yml"      ]]; then TEMPLATE_YML="$RUN_DIR/taiji_config.template.yml"
else echo "ERROR: taiji_config.template.yml missing (looked in $RUN_DIR/code/ and $RUN_DIR/)" >&2; exit 2
fi

echo "[run-taiji] REPO_ROOT=$REPO_ROOT"
echo "[run-taiji] RUN_DIR=$RUN_DIR"
echo "[run-taiji] samples=$SAMPLES_CSV"
echo "[run-taiji] template=$TEMPLATE_YML"
echo "[run-taiji] runner args: ${EXTRA_RUNNER_ARGS[*]}"

# ---- Step 1: pick the binary --------------------------------------------

if [[ -n "$BINARY_OVERRIDE" ]]; then
  BINARY="$BINARY_OVERRIDE"
else
  case "$(uname -s)-$(uname -m)" in
    Darwin-*)
      BINARY="$REPO_ROOT/binaries/taiji-macOS-Catalina-10.15"
      ;;
    Linux-x86_64)
      if [[ -f "$REPO_ROOT/binaries/taiji-Ubuntu-x86_64" ]]; then
        BINARY="$REPO_ROOT/binaries/taiji-Ubuntu-x86_64"
      else
        BINARY="$REPO_ROOT/binaries/taiji-CentOS-x86_64"
      fi
      ;;
    *)
      if [[ "$DRY_RUN" -eq 1 ]] || [[ "$PREPARE_ONLY" -eq 1 ]]; then
        # Under --dry-run / --prepare-only we don't actually invoke the binary;
        # pick any existing one as a placeholder for the planned command.
        BINARY=$(ls -1 "$REPO_ROOT/binaries/" 2>/dev/null | head -1)
        BINARY="$REPO_ROOT/binaries/${BINARY:-none}"
        echo "[run-taiji] WARNING: $(uname -sm) has no native Taiji binary;" >&2
        echo "[run-taiji] using $BINARY as a placeholder (no-op for this run)." >&2
      else
        echo "ERROR: unsupported OS+arch '$(uname -sm)'." >&2
        echo "Pass --binary <path> to override." >&2
        exit 3
      fi
      ;;
  esac
fi

[[ -f "$BINARY" ]] || { echo "ERROR: binary not found: $BINARY" >&2; exit 3; }
chmod +x "$BINARY" 2>/dev/null || true
echo "[run-taiji] binary=$BINARY"

PYTHON="${PYTHON:-python3}"
command -v "$PYTHON" >/dev/null 2>&1 || { echo "ERROR: $PYTHON not on PATH" >&2; exit 3; }

# ---- Step 2: regenerate xlsx with absolute paths ------------------------

echo
echo "[run-taiji] regenerating $RUN_DIR/taiji_input.xlsx ..."
"$PYTHON" "$REPO_ROOT/skills/build-taiji-input/scripts/build_taiji_input.py" \
  --samples  "$SAMPLES_CSV" \
  --data-dir "$REPO_ROOT/data/data" \
  --genome   hg38 \
  --rna-pattern  '{group}_RNA.tsv' \
  --atac-pattern '{group}_REP1.mLb.clN_peaks_small.narrowPeak' \
  --hic-pattern  '{group}_hicloops_chr22.bedpe' \
  --out      "$RUN_DIR/taiji_input.xlsx"

# ---- Step 3: pre-flight -------------------------------------------------

if [[ "$SKIP_PREFLIGHT" -eq 0 ]]; then
  echo
  echo "[run-taiji] running pre-flight ..."
  if ! "$PYTHON" "$REPO_ROOT/bin/preflight-xlsx.py" "$RUN_DIR/taiji_input.xlsx"; then
    echo "[run-taiji] PREFLIGHT FAILED — refusing to launch Taiji." >&2
    echo "[run-taiji] Pass --skip-preflight to override (not recommended)." >&2
    exit 4
  fi
fi

# ---- Step 4: start workflow-log -----------------------------------------

LOG_PY="$REPO_ROOT/skills/workflow-log/scripts/log.py"
REF_JSON=$(printf '{"fasta":"%s/dependencies_data/hg38/genome.fa","gtf":"%s/dependencies_data/hg38/genes.gtf","motif":"%s/cisbp_database/cisbp_human_2.meme"}' \
  "$REPO_ROOT" "$REPO_ROOT" "$REPO_ROOT")

echo
echo "[run-taiji] starting workflow-log ..."
RUN_ID=$("$PYTHON" "$LOG_PY" start \
  --work-dir "$RUN_DIR" \
  --genome hg38 \
  --reference-paths "$REF_JSON" \
  --note "via bin/run-taiji.sh on $(uname -sm)")
echo "[run-taiji] RUN_ID=$RUN_ID"

FINALIZED=0
WALL_DURATION=0
finalize_log() {
  local status="$1"
  local note="${2:-}"
  if [[ "$FINALIZED" -eq 0 ]]; then
    "$PYTHON" "$LOG_PY" finalize \
      --work-dir "$RUN_DIR" --status "$status" \
      ${note:+--note "$note"} \
      --duration-s "$WALL_DURATION" || true
    FINALIZED=1
  fi
}
trap 'finalize_log failure "wrapper exited unexpectedly"' EXIT
trap 'finalize_log failure "interrupted"; exit 130' INT TERM

# ---- Step 5: delegate to taiji-runner skill -----------------------------

RUNNER_PY="$REPO_ROOT/skills/taiji-runner/scripts/run_taiji.py"

RUNNER_ARGS=(
  --xlsx     "$RUN_DIR/taiji_input.xlsx"
  --template "$TEMPLATE_YML"
  --run-dir  "$RUN_DIR"
  --binary   "$BINARY"
  --repo-root "$REPO_ROOT"
)
[[ -n "$SAMPLES" ]]            && RUNNER_ARGS+=(--samples "$SAMPLES")
[[ "$PREPARE_ONLY" -eq 1 ]]    && RUNNER_ARGS+=(--prepare-only)
[[ "$DRY_RUN"      -eq 1 ]]    && RUNNER_ARGS+=(--prepare-only)   # dry-run -> prepare-only
RUNNER_ARGS+=("${EXTRA_RUNNER_ARGS[@]}")

# Wrap the runner in sandbox-run.sh so all Taiji writes stay in the workspace
# and TMPDIR points at REPO_ROOT/tmp (instead of /tmp). --no-sandbox bypasses
# the wrapper entirely (e.g. for debugging). --strict-sandbox refuses to run
# unless an OS-level sandbox (bwrap on Linux, sandbox-exec on macOS) is
# available.
SANDBOX_PREFIX=()
if [[ "$NO_SANDBOX" -eq 0 ]]; then
  SANDBOX_PREFIX=("$REPO_ROOT/bin/sandbox-run.sh")
  [[ "$STRICT_SANDBOX" -eq 1 ]] && SANDBOX_PREFIX+=(--strict)
  SANDBOX_PREFIX+=(--)
fi

echo
echo "[run-taiji] command: ${SANDBOX_PREFIX[*]:-} $PYTHON $RUNNER_PY ${RUNNER_ARGS[*]}"
START_TS=$SECONDS
RUNNER_EXIT=0
"${SANDBOX_PREFIX[@]}" "$PYTHON" "$RUNNER_PY" "${RUNNER_ARGS[@]}" || RUNNER_EXIT=$?
WALL_DURATION=$((SECONDS - START_TS))

echo
echo "[run-taiji] taiji-runner exited $RUNNER_EXIT  duration=${WALL_DURATION}s"

# ---- Step 6: finalize ---------------------------------------------------

if [[ "$RUNNER_EXIT" -eq 0 ]]; then
  finalize_log success "all samples ok"
else
  finalize_log failure "taiji-runner exited $RUNNER_EXIT — see log/ for per-sample details"
fi
trap - EXIT INT TERM

if [[ "$DRY_RUN" -eq 1 || "$PREPARE_ONLY" -eq 1 ]]; then
  echo "[run-taiji] DONE (prepare/dry-run mode). Inspect $RUN_DIR/Input/ + Output/ tree."
  exit 0
fi

if [[ "$RUNNER_EXIT" -eq 0 ]]; then
  echo "[run-taiji] DONE. Outputs: $RUN_DIR/Output/Partial/"
  echo "[run-taiji] Log: $RUN_DIR/log/${RUN_ID}.md"
  exit 0
else
  echo "[run-taiji] FAILED. Inspect per-sample stderr in $RUN_DIR/log/<sample>.taiji.stderr"
  exit 5
fi
