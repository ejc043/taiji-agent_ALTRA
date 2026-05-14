#!/bin/bash
#SBATCH --job-name=d7_chooseR
#SBATCH --output=/stg3/data1/eunice/Projects/taiji_agent/small_test/taiji-agent/runs/d7-siiel-spleen-mm10/coembed/chooseR_%j.log
#SBATCH --time=4:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu
#SBATCH --dependency=

REPO=/stg3/data1/eunice/Projects/taiji_agent/small_test/taiji-agent
RSCRIPT=/stg3/data1/eunice/.local/share/mamba/envs/taiji-agent/bin/Rscript

cd "$REPO"

$RSCRIPT skills/coembed-construct/scripts/chooseR_clustering.R \
  --input              runs/d7-siiel-spleen-mm10/coembed/coembed.rds \
  --output             runs/d7-siiel-spleen-mm10/coembed/coembed_clustered.rds \
  --resolutions        0.1,0.2,0.4,0.6,0.8,1.0,1.2,1.5,2.0,2.5,3.0,4.0,5.0 \
  --n-subsample        2000 \
  --n-bootstraps       10 \
  --n-pcs              30 \
  --min-cluster-cells  20 \
  --metadata-cols      tissue,genotype \
  --seed               42
