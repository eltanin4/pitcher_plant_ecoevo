#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 6:00:00
#SBATCH --array=9-10
#SBATCH -p newnodes,sched_mit_chisholm,sched_mit_hill,sched_any
#SBATCH -J qntbin
#SBATCH --mem=50G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/refinement_%a.out

# # Setting checkm data directory.
# checkm data setRoot /home/akshitg/checkm_data_dir

# Making the bin refinement directory.
rm -r /home/akshitg/corderolab/akshitg/Binning/M0"${SLURM_ARRAY_TASK_ID}"-binning/QUANT_BINS/
mkdir /home/akshitg/corderolab/akshitg/Binning/M0"${SLURM_ARRAY_TASK_ID}"-binning/QUANT_BINS

# Quantifying bin abundances for identified bins.
metaWRAP quant_bins \
    -b /home/akshitg/corderolab/akshitg/Binning/M0"${SLURM_ARRAY_TASK_ID}"-binning/BIN_REFINEMENT/metawrap_50_10_bins \
    -o /home/akshitg/corderolab/akshitg/Binning/M0"${SLURM_ARRAY_TASK_ID}"-binning/QUANT_BINS \
    -a /home/akshitg/corderolab/akshitg/Assemblies/M0"${SLURM_ARRAY_TASK_ID}"-assembly/out/final.contigs.fa \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260*.fastq
