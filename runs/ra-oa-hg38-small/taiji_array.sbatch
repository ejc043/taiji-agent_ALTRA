#!/bin/bash
# taiji_array.sbatch — SLURM array script for per-sample Taiji execution.
#
# This file is a TEMPLATE materialized at submit time by the taiji-runner
# skill (skills/taiji-runner/scripts/run_taiji.py). It substitutes six
# tokens for: job name, log dir, cpus per task, memory in GB, wall-time
# limit, and the absolute path to the Taiji binary. See the SKILL.md /
# run_taiji.py source for the exact substitution map.
#
# Invoked by sbatch as:
#   sbatch -a 1-N%MAX <this-script> <path-to-taiji_config_files.txt>
#
# Each array task picks the Nth config from the list, cd's into the
# matching per-sample output dir, and runs `taiji run --config <yml>`.

#SBATCH --job-name=taiji-ra-oa-hg38-small
#SBATCH --output=/stg3/data1/eunice/Projects/taiji_agent/small_test/taiji-agent/runs/ra-oa-hg38-small/log/slurm-%A_%a.out
#SBATCH --error=/stg3/data1/eunice/Projects/taiji_agent/small_test/taiji-agent/runs/ra-oa-hg38-small/log/slurm-%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=04:00:00

set -euo pipefail

CONFIG_LIST="${1:?usage: sbatch <script> <taiji_config_files.txt>}"
CONFIG=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CONFIG_LIST")
[[ -z "$CONFIG" ]] && { echo "no config for array task $SLURM_ARRAY_TASK_ID" >&2; exit 2; }
[[ -f "$CONFIG" ]] || { echo "config not found: $CONFIG" >&2; exit 2; }

# Resolve the per-sample output_dir from the config (Taiji writes there).
# `output_dir:` is one of the config's first fields; awk picks it out.
OUTPUT_DIR=$(awk -F': *' '/^output_dir:/ {print $2; exit}' "$CONFIG" | tr -d '"' | tr -d "'")
[[ -d "$OUTPUT_DIR" ]] || mkdir -p "$OUTPUT_DIR"

echo "[task ${SLURM_ARRAY_TASK_ID}] config=$CONFIG"
echo "[task ${SLURM_ARRAY_TASK_ID}] cwd=$OUTPUT_DIR"
echo "[task ${SLURM_ARRAY_TASK_ID}] cpus=${SLURM_CPUS_PER_TASK}  mem=${SLURM_MEM_PER_NODE:-?}"
echo "[task ${SLURM_ARRAY_TASK_ID}] binary=/stg3/data1/eunice/Projects/taiji_agent/small_test/taiji-agent/binaries/taiji-CentOS-x86_64"

cd "$OUTPUT_DIR"

# +RTS -N<n>: Taiji is Haskell; thread count is an RTS option, not a CLI flag.
/stg3/data1/eunice/Projects/taiji_agent/small_test/taiji-agent/binaries/taiji-CentOS-x86_64 run --config "$CONFIG" +RTS "-N${SLURM_CPUS_PER_TASK}"
