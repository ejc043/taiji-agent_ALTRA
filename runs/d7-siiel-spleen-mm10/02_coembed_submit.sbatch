#!/bin/bash
#SBATCH --job-name=d7_coembed
#SBATCH --output=/stg3/data1/eunice/Projects/taiji_agent/small_test/taiji-agent/runs/d7-siiel-spleen-mm10/coembed/coembed_%j.log
#SBATCH --time=4:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu
#SBATCH --dependency=afterok:872786

REPO=/stg3/data1/eunice/Projects/taiji_agent/small_test/taiji-agent
RSCRIPT=/stg3/data1/eunice/.local/share/mamba/envs/taiji-agent/bin/Rscript

cd "$REPO"

mkdir -p runs/d7-siiel-spleen-mm10/coembed

$RSCRIPT skills/coembed-construct/scripts/coembed.R \
  --rna           runs/d7-siiel-spleen-mm10/qc/rna_filtered.rds \
  --atac          data/singlecell/D7_scATAC_prep.rds \
  --fragments     data/singlecell/D7_fragments.tsv.gz \
  --genome        mm10 \
  --output        runs/d7-siiel-spleen-mm10/coembed/coembed.rds \
  --metadata-cols tissue,genotype \
  --reuse-atac-reductions \
  --n-pcs         30 \
  --no-cluster
