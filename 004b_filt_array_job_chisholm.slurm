#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --array=46-47
#SBATCH -t 20:00:00
#SBATCH -p sched_mit_chisholm
#SBATCH -J 10_filt_uniq
#SBATCH --mem=50G

~/.linuxbrew/bin/python3.7 004a1_filter_reads_mapq_60_longest.py $SLURM_ARRAY_TASK_ID
