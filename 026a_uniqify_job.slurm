#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 36:00:00
#SBATCH -p sched_mit_chisholm
#SBATCH -J unq_fb
#SBATCH --mem=100G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/unq_snps_freebayes_no_threshold.out

python3 026_uniqify_snps.py
