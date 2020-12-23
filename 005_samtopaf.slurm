#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --array=100-114
#SBATCH -t 12:00:00
#SBATCH -p newnodes,sched_mit_chisholm,sched_mit_hill,sched_any
#SBATCH -J sam_paf_100
#SBATCH --mem=20G

~/minimap2-2.17_x64-linux/paftools.js sam2paf "D19-260"$SLURM_ARRAY_TASK_ID"_minimap2_axsr_results.sam" > "D19-260"$SLURM_ARRAY_TASK_ID"_minimap2_axsr_results.paf"
