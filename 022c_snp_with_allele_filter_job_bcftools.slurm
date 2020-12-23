#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 18:00:00
#SBATCH --array=35-114
#SBATCH -p sched_mit_chisholm
#SBATCH -J proc_snp_bcf
#SBATCH --mem=60G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/proc_snp_bcftools_%a.out

for ((i=1; i<=24; i++));
    do
        python3 /home/akshitg/lte/metagenomes_lte/scripts/022c_filter_snps_with_allele_freqs_bcftools.py $SLURM_ARRAY_TASK_ID $i
    done
for ((i=26; i<=34; i++));
    do
        python3 /home/akshitg/lte/metagenomes_lte/scripts/022c_filter_snps_with_allele_freqs_bcftools.py $SLURM_ARRAY_TASK_ID $i
    done
