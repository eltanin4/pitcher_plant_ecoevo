#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --array=35-114
#SBATCH -t 6:00:00
#SBATCH -p newnodes,sched_mit_chisholm,sched_mit_hill,sched_any
#SBATCH -J avg_cov
#SBATCH --mem=64G

~/.linuxbrew/bin/python3.7 ~/lte/metagenomes_lte/scripts/010_avg_cov_calc.py $SLURM_ARRAY_TASK_ID
