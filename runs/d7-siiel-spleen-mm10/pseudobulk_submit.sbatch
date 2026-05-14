#!/bin/bash
#SBATCH --job-name=d7_pseudobulk
#SBATCH --output=/stg3/data1/eunice/Projects/taiji_agent/small_test/taiji-agent/runs/d7-siiel-spleen-mm10/pseudobulk/pseudobulk_%j.log
#SBATCH --time=8:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --partition=cpu
#SBATCH --dependency=afterok:872422

REPO=/stg3/data1/eunice/Projects/taiji_agent/small_test/taiji-agent
PYTHON=/stg3/data1/eunice/.local/share/mamba/envs/taiji-agent/bin/python

cd "$REPO"

mkdir -p runs/d7-siiel-spleen-mm10/pseudobulk

$PYTHON skills/pseudobulk-construct/scripts/pseudobulk.py \
  --input         runs/d7-siiel-spleen-mm10/coembed/coembed_clustered.rds \
  --fragments     data/singlecell/D7_fragments.tsv.gz \
  --output-dir    runs/d7-siiel-spleen-mm10/pseudobulk \
  --metadata-cols tissue,genotype \
  --cohort-col    tissue \
  --use-existing-clusters \
  --genome        mm10
