#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH -p sched_mit_chisholm
#SBATCH -J unq_snp
#SBATCH --mem=100G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/unq_snps_bcftools.out

python3 028d_uniqify_snps.py
