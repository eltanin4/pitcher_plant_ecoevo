#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 18:00:00
#SBATCH --array=1-10
#SBATCH -p sched_mit_chisholm
#SBATCH -J traj
#SBATCH --mem=50G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/trajectories_no_threshold_%a.out

python3 033c_snp_trajectories_no_threshold.py $SLURM_ARRAY_TASK_ID
