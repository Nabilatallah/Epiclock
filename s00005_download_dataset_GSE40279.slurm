#!/bin/bash
#SBATCH --job-name=s00005_download_dataset_GSE40279
#SBATCH --output=/home/n.atallah/epiclock/ec01/logs/s00005_download_dataset_GSE40279_%j.out
#SBATCH --error=/home/n.atallah/epiclock/ec01/logs/s00005_download_dataset_GSE40279_%j.err
#SBATCH --time=04:00:00
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=n.atallah@northeastern.edu

echo "🚀 Starting job: download_gse40279 at $(date)"

# === Run Python script inside container with conda env activated ===
apptainer exec \
    --bind /home/n.atallah/epiclock:/mnt \
    /home/n.atallah/epiclock/containers/epiclock_env.sif \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
             conda activate ml_em_env && \
             python /home/n.atallah/epiclock/ec01/code/s00005_download_dataset_GSE40279.py \
             --geo_id GSE40279 --sample_index 0"

echo "✅ Job completed at $(date)"
