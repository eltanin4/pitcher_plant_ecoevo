#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --array=39-40
#SBATCH -t 20:00:00
#SBATCH -p newnodes,sched_mit_chisholm,sched_mit_hill,sched_any
#SBATCH -J 10_filt_uniq
#SBATCH --mem=60G

~/.linuxbrew/bin/python3.7 004a1_filter_reads_mapq_60_longest.py $SLURM_ARRAY_TASK_ID
