#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 10
#SBATCH -t 16:00:00
#SBATCH -p sched_mit_chisholm
#SBATCH -J dNdS
#SBATCH --mem=60G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/dNdS.out

python3 027_calculate_dnds_genome_wise.py
