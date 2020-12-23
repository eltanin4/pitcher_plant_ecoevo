#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH --array=1-10
#SBATCH -p newnodes,sched_mit_chisholm,sched_mit_hill,sched_any
#SBATCH -J traj
#SBATCH --mem=50G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/trajectories_%a.out

python3 033_snp_trajectories_with_threshold_bt2.py $SLURM_ARRAY_TASK_ID

