#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH --array=1-34
#SBATCH -p newnodes,sched_mit_chisholm,sched_any,sched_mit_hill
#SBATCH -J snp_corr
#SBATCH --mem=30G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/snp_correlations_spearman_%a.out

python3 036_snp_correlations_spearman_bt2.py $SLURM_ARRAY_TASK_ID
