#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --array=35-99
#SBATCH -t 1:00:00
#SBATCH -p newnodes,sched_mit_chisholm,sched_mit_hill,sched_any
#SBATCH -J cov10
#SBATCH --mem=30G

bedtools genomecov -ibam ~/lte/metagenomes_lte/D19-2600${SLURM_ARRAY_TASK_ID}_minimap2_axsr_results.bam -bg > ~/lte/metagenomes_lte/bgraph_cov_mg_0${SLURM_ARRAY_TASK_ID}.tsv
