#!/bin/bash
#SBATCH --array=1-14

source ~/.bashrc
conda activate imlabtools
echo "Job starts at $(date)"

mkdir -p ~/airr/results/assoc_pheno/$2
bash ~/airr/scripts/ukb/one_assoc.sh ${SLURM_ARRAY_TASK_ID} $1 $2

echo "Job ends at $(date)"