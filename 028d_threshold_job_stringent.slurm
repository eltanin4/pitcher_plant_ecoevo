#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 18:00:00
#SBATCH --array=1-10
#SBATCH -p sched_mit_chisholm
#SBATCH -J threshold
#SBATCH --mem=50G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/threshold_%a.out

python3 028d_threshold_snps_per_microcosm_stringent.py $SLURM_ARRAY_TASK_ID
