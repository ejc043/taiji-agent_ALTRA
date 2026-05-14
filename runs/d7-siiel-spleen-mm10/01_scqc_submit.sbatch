#!/bin/bash
#SBATCH --job-name=d7_scqc
#SBATCH --output=/stg3/data1/eunice/Projects/taiji_agent/small_test/taiji-agent/runs/d7-siiel-spleen-mm10/qc/scqc_%j.log
#SBATCH --time=0:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --partition=cpu

REPO=/stg3/data1/eunice/Projects/taiji_agent/small_test/taiji-agent
RSCRIPT=/stg3/data1/eunice/.local/share/mamba/envs/taiji-agent/bin/Rscript

cd "$REPO"

mkdir -p runs/d7-siiel-spleen-mm10/qc

$RSCRIPT skills/sc-qc/scripts/filter_rna.R \
  --input        data/singlecell/D7_scRNA_prep.rds \
  --output       runs/d7-siiel-spleen-mm10/qc/rna_filtered.rds \
  --min-features 200 \
  --max-features 5000 \
  --max-pct-mt   10 \
  --mt-pattern   "^mt-" \
  --group-by     tissue
