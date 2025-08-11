#!/bin/bash
#SBATCH --job-name=s00007_extract_annotation.py
#SBATCH --output=/home/n.atallah/epiclock/ec01/logs/s00007_extract_annotation_%j.out
#SBATCH --error=/home/n.atallah/epiclock/ec01/logs/s00007_extract_annotation_%j.err
#SBATCH --time=04:00:00
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=n.atallah@northeastern.edu

echo "🚀 Starting s00007_extract_annotation.py at $(date)"

# === Run Python script inside container with conda env activated ===
apptainer exec \
    --bind /home/n.atallah/epiclock:/mnt \
    /home/n.atallah/epiclock/containers/epiclock_env.sif \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate ml_em_env && \
             python /mnt/ec01/code/s00007_extract_annotation.py \
               --process-rda \
               --tag Illumina450k \
               --data_dir /mnt/ec01/data \
               --results_dir /mnt/ec01/results \
               --log-level INFO"

echo "✅ Job completed at $(date)"
