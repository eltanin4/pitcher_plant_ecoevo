#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH --array=35-114
#SBATCH -p newnodes,sched_mit_chisholm,sched_mit_hill
#SBATCH -J gene_snp
#SBATCH --mem=50G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/gene_wise_snps_%a.out

python3 /home/akshitg/lte/metagenomes_lte/scripts/031c_per_time_gene_wise_snp_table.py $SLURM_ARRAY_TASK_ID
