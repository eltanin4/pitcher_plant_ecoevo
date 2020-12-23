#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 20
#SBATCH -t 72:00:00
#SBATCH --array=5-6
#SBATCH -p sched_mit_chisholm
#SBATCH -J reassemble
#SBATCH --mem=100G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/reassemble_%a.out

# # Setting checkm data directory.
# checkm data setRoot /home/akshitg/checkm_data_dir

# Making the bin refinement directory.
rm -r /home/akshitg/corderolab/akshitg/Binning/M0"${SLURM_ARRAY_TASK_ID}"-binning/BIN_REASSEMBLY/
mkdir /home/akshitg/corderolab/akshitg/Binning/M0"${SLURM_ARRAY_TASK_ID}"-binning/BIN_REASSEMBLY

# Concatenating all relevant microcosm reads into two giant fastq files.
echo "Concatenating all reads"
cat /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-2600*_1.fastq \
    > /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_ALL_READS_1.fastq
cat /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-2600*_2.fastq \
    > /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_ALL_READS_2.fastq

# Re-assembling bins using the refined bins.
echo "Re-assembling bins"
metaWRAP reassemble_bins \
    -o /home/akshitg/corderolab/akshitg/Binning/M0"${SLURM_ARRAY_TASK_ID}"-binning/BIN_REASSEMBLY \
    -1 /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_ALL_READS_1.fastq \
    -2 /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_ALL_READS_2.fastq \
    -t 20 \
    -m 96 \
    -c 50 \
    -x 10 \
    -b /home/akshitg/corderolab/akshitg/Binning/M0"${SLURM_ARRAY_TASK_ID}"-binning/BIN_REFINEMENT/metawrap_50_10_bins

# Removing files
echo "Removing files"
rm /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_ALL_READS_*.fastq
