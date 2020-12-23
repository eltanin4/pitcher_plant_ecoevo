#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 20
#SBATCH -t 12:00:00
#SBATCH --array=1-7
#SBATCH -p newnodes,sched_mit_chisholm,sched_mit_hill,sched_any
#SBATCH -J refine
#SBATCH --mem=100G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/refinement_%a.out

# Setting checkm data directory.
checkm data setRoot /home/akshitg/checkm_data_dir

# Making the bin refinement directory.
rm -r /home/akshitg/corderolab/akshitg/Binning/M0"${SLURM_ARRAY_TASK_ID}"-binning/BIN_REFINEMENT/
mkdir /home/akshitg/corderolab/akshitg/Binning/M0"${SLURM_ARRAY_TASK_ID}"-binning/BIN_REFINEMENT

# Running bin refinment with a completion threshold of >= 50% and contamination threshold of <= 10%.
metawrap bin_refinement \
    -o /home/akshitg/corderolab/akshitg/Binning/M0"${SLURM_ARRAY_TASK_ID}"-binning/BIN_REFINEMENT \
    -t 20 \
    -m 100 \
    -A /home/akshitg/corderolab/akshitg/Binning/M0"${SLURM_ARRAY_TASK_ID}"-binning/INITIAL_BINNING/metabat2_bins/ \
    -B /home/akshitg/corderolab/akshitg/Binning/M0"${SLURM_ARRAY_TASK_ID}"-binning/INITIAL_BINNING/maxbin2_bins/ \
    -C /home/akshitg/corderolab/akshitg/Binning/M0"${SLURM_ARRAY_TASK_ID}"-binning/INITIAL_BINNING/concoct_bins/ \
    -c 50 -x 10
