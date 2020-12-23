#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH --array=5
#SBATCH -p sched_mit_chisholm
#SBATCH -J traj
#SBATCH --mem=50G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/trajectories_%a.out

python3 033_snp_trajectories_from_processed_snps.py $SLURM_ARRAY_TASK_ID
