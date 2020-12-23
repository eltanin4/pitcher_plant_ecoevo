#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 20
#SBATCH -t 18:00:00
#SBATCH -p sched_mit_chisholm
#SBATCH -J gene_snp
#SBATCH --mem=100G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/gene_wise_snps.out

python3 031_gene_wise_snp_table.py
