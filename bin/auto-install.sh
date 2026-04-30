#!/usr/bin/env bash
# auto-install.sh — agent-driven install for the taiji-agent.
#
# Runs `detect-dataset-type` on the supplied --data-dir, picks the right
# profile based on the classification, and execs `bin/install.sh` with that
# profile. Lets users get a working install with one command instead of
# remembering "bulk -> --profile base, single-cell -> --profile sc".
#
# Decision rules:
#   classification == "bulk"        -> profile=base
#   classification == "single-cell" -> profile=sc
#   classification == "mixed"       -> profile=sc (since either branch may be used)
#   classification == "unknown"     -> profile=base + warning (fall back to safe minimum)
#
# Usage:
#   bash bin/auto-install.sh --data-dir <path> [--system macos|centos|ubuntu]
#                            [--solver micromamba|mamba|conda]
#                            [--fetch-references --genome hg38|hg19|mm10]
#                                                       # stage reference data in the same step
#                            [--dry-run]                # print decision, don't install
#                            [...any other install.sh flag is forwarded]
#
# Examples:
#   bash bin/auto-install.sh --data-dir data/RA_OA_chr22_demo --system macos
#   bash bin/auto-install.sh --data-dir data/RA_OA_chr22_demo --system macos \
#                            --fetch-references --genome hg38
#   bash bin/auto-install.sh --data-dir data/some_pbmc_multiome --system centos --dry-run

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DETECT_PY="$REPO_ROOT/skills/detect-dataset-type/scripts/detect.py"
INSTALL_SH="$REPO_ROOT/bin/install.sh"

DATA_DIR=""
DRY_RUN=0
PASS_THROUGH=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-dir) DATA_DIR="$2"; shift 2 ;;
    --dry-run)  DRY_RUN=1; shift ;;
    -h|--help)  sed -n '2,25p' "$0"; exit 0 ;;
    *)          PASS_THROUGH+=("$1"); shift ;;
  esac
done

if [[ -z "$DATA_DIR" ]]; then
  echo "ERROR: --data-dir <path> is required (the dataset to classify)" >&2
  echo "       e.g. bash bin/auto-install.sh --data-dir data/my_dataset --system macos" >&2
  exit 2
fi
if [[ ! -d "$DATA_DIR" ]]; then
  echo "ERROR: --data-dir $DATA_DIR does not exist or is not a directory" >&2
  exit 2
fi

PYTHON="${PYTHON:-python3}"
command -v "$PYTHON" >/dev/null 2>&1 || { echo "ERROR: $PYTHON not on PATH" >&2; exit 3; }

# ---- Run detect-dataset-type ----
# Capture stderr to a temp file so any error from detect.py is visible if
# stdout comes back empty (otherwise `2>/dev/null` would hide the real cause).
echo "[auto-install] classifying $DATA_DIR ..."
DETECT_STDERR=$(mktemp)
DETECT_JSON=$("$PYTHON" "$DETECT_PY" "$DATA_DIR" --format json 2>"$DETECT_STDERR" || true)
if [[ -z "$DETECT_JSON" ]]; then
  echo "ERROR: detect-dataset-type returned no output" >&2
  if [[ -s "$DETECT_STDERR" ]]; then
    echo "  stderr from detect:" >&2
    sed 's/^/    /' "$DETECT_STDERR" >&2
  else
    echo "  (no stderr either — try running detect.py manually:" >&2
    echo "     $PYTHON $DETECT_PY $DATA_DIR --format json)" >&2
  fi
  rm -f "$DETECT_STDERR"
  exit 3
fi
rm -f "$DETECT_STDERR"

CLASSIFICATION=$("$PYTHON" -c "import json,sys; print(json.loads(sys.stdin.read()).get('classification','unknown'))" <<<"$DETECT_JSON")
SC_MODALITY=$("$PYTHON"   -c "import json,sys; print(json.loads(sys.stdin.read()).get('sc_modality') or '')"      <<<"$DETECT_JSON")

echo "[auto-install] classification: $CLASSIFICATION"
[[ -n "$SC_MODALITY" ]] && echo "[auto-install] sc_modality:    $SC_MODALITY"

# ---- Pick profile ----
case "$CLASSIFICATION" in
  bulk)
    PROFILE="base"
    REASON="bulk dataset — only base profile (Python + xlsx + MACS3) needed"
    ;;
  single-cell)
    PROFILE="sc"
    REASON="single-cell dataset — base + sc (R + Seurat + Signac) needed"
    ;;
  mixed)
    PROFILE="sc"
    REASON="mixed dataset — installing base + sc to cover both branches"
    ;;
  unknown)
    PROFILE="base"
    REASON="WARNING: classification was 'unknown'; falling back to base. "\
"If you actually have SC data, re-run with: bash bin/install.sh --profile sc"
    ;;
  *)
    echo "ERROR: unexpected classification '$CLASSIFICATION'" >&2
    exit 3
    ;;
esac

echo "[auto-install] decision: --profile $PROFILE"
echo "[auto-install] reason:   $REASON"
echo

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "[auto-install] --dry-run set; not invoking install.sh"
  echo "[auto-install] would run: bash bin/install.sh --profile $PROFILE ${PASS_THROUGH[*]:-}"
  exit 0
fi

# ---- Forward to install.sh ----
echo "[auto-install] running: bash bin/install.sh --profile $PROFILE ${PASS_THROUGH[*]:-}"
exec bash "$INSTALL_SH" --profile "$PROFILE" "${PASS_THROUGH[@]}"
