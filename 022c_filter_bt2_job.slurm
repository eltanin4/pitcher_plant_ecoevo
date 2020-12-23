#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH --array=35-114
#SBATCH -p newnodes,sched_mit_chisholm,sched_any,sched_mit_hill
#SBATCH -J proc_snp
#SBATCH --mem=30G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/proc_snp_%a.out

for ((i=1; i<=24; i++));
    do
        python3 /home/akshitg/lte/metagenomes_lte/scripts/022c_filter_snps_with_allele_freqs_bt2.py $SLURM_ARRAY_TASK_ID $i
    done
for ((i=26; i<=34; i++));
    do
        python3 /home/akshitg/lte/metagenomes_lte/scripts/022c_filter_snps_with_allele_freqs_bt2.py $SLURM_ARRAY_TASK_ID $i
    done
