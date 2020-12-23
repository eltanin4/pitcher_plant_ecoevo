#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH --array=1-10
#SBATCH -p newnodes,sched_mit_chisholm,sched_mit_hill,sched_any
#SBATCH -J threshold
#SBATCH --mem=50G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/bcftools_threshold_%a.out

python3 028c_threshold_snps_per_microcosm.py $SLURM_ARRAY_TASK_ID
