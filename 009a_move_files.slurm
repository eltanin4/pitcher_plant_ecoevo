#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --array=35-99
#SBATCH -t 3:00:00
#SBATCH -p newnodes,sched_mit_chisholm,sched_mit_hill,sched_any
#SBATCH -J mover
#SBATCH --mem=100G

mv -v ~/lte/metagenomes_lte/D19-2600${SLURM_ARRAY_TASK_ID}_minimap2_axsr_results.sam ~/main/lte/minimap2_sams
