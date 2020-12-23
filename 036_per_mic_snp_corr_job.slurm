#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 18:00:00
#SBATCH --array=1-34
#SBATCH -p sched_mit_chisholm
#SBATCH -J snp_corr
#SBATCH --mem=50G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/snp_correlations_%a.out

python3 036_snp_correlations.py $SLURM_ARRAY_TASK_ID
