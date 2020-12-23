#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 10:00:00
#SBATCH --array=74
#SBATCH -p newnodes,sched_mit_chisholm,sched_mit_hill,sched_any
#SBATCH -J proc_snp
#SBATCH --mem=50G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/proc_snp_%a.out

for ((i=1; i<=24; i++));
    do
        python3 /home/akshitg/lte/metagenomes_lte/scripts/022_filter_snps_from_tsvs_freebayes.py $SLURM_ARRAY_TASK_ID $i
    done
for ((i=26; i<=34; i++));
    do
        python3 /home/akshitg/lte/metagenomes_lte/scripts/022_filter_snps_from_tsvs_freebayes.py $SLURM_ARRAY_TASK_ID $i
    done
