#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --array=39
#SBATCH -t 2:00:00
#SBATCH -p newnodes,sched_mit_chisholm,sched_mit_hill,sched_any
#SBATCH -J uniq_reads
#SBATCH --mem=40G

~/.linuxbrew/bin/python3.7 006a1_uniqify_reads.py $SLURM_ARRAY_TASK_ID


