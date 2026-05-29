#!/bin/bash
# bin/setup_altra.sh
# Creates the taiji-agent-altra conda environment and installs ArchR.
# Called by Claude when the user asks to set up or make the ALTRA environment.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_NAME="taiji-agent-altra"
ENV_FILE="${REPO_ROOT}/environment.altra.yml"
POSTINSTALL="${REPO_ROOT}/bin/postinstall_altra.R"

# Detect micromamba or conda
if command -v micromamba &>/dev/null; then
  CONDA_CMD="micromamba"
elif command -v conda &>/dev/null; then
  CONDA_CMD="conda"
else
  echo "ERROR: neither micromamba nor conda found in PATH." >&2
  echo "Install micromamba: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html" >&2
  exit 1
fi

echo "Using: $($CONDA_CMD --version 2>&1 | head -1)"
echo ""

# Check if env already exists
if $CONDA_CMD env list 2>/dev/null | grep -q "^${ENV_NAME}[[:space:]]"; then
  echo "Environment '${ENV_NAME}' already exists."
  SKIP_CREATE=1
else
  SKIP_CREATE=0
fi

# Step 1: create env
if [ "$SKIP_CREATE" -eq 0 ]; then
  echo "=== Step 1/2: Creating ${ENV_NAME} from environment.altra.yml ==="
  echo "(This takes ~20-40 min on first run)"
  $CONDA_CMD create -f "${ENV_FILE}" -y
  echo "Environment created."
else
  echo "=== Step 1/2: Skipping env creation (already exists) ==="
fi

# Step 2: install ArchR via postinstall_altra.R
echo ""
echo "=== Step 2/2: Installing ArchR from GitHub ==="
RSCRIPT=$($CONDA_CMD run -n "${ENV_NAME}" which Rscript 2>/dev/null)
if [ -z "$RSCRIPT" ]; then
  echo "ERROR: Rscript not found in ${ENV_NAME}." >&2
  exit 1
fi

$CONDA_CMD run -n "${ENV_NAME}" Rscript "${POSTINSTALL}"

echo ""
echo "=== Setup complete ==="
echo "Activate with:  micromamba activate ${ENV_NAME}"
echo "Then run:       claude"
