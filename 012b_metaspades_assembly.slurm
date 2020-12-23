#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --array=9-10
#SBATCH -t 96:00:00
#SBATCH -p sched_mit_chisholm
#SBATCH -J metaspades
#SBATCH --mem=100G

# SLURM_ARRAY_TASK_ID corresponds to microcosm number.
# NOTE: M09 and M10 have to be dealt with differently 
# because of their numbering schemes.
# SLURM_ARRAY_TASK_ID=1

# Adding symbolic links to all fastq files.
echo "Creating symbolic links"
for ((i=0; i<8; i++)); do
    ln -s \
    /home/akshitg/corderolab/akshitg/D19-2600"$((90 + SLURM_ARRAY_TASK_ID + i * 2))"-4151H/190801Cor_D19-2600"$(( \
    90 + SLURM_ARRAY_TASK_ID + i * 2))"_1_sequence.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-2600"$(( \
    90 + SLURM_ARRAY_TASK_ID + i * 2))"_1.fastq
    ln -s \
    /home/akshitg/corderolab/akshitg/D19-2600"$((90 + SLURM_ARRAY_TASK_ID + i * 2))"-4151H/190801Cor_D19-2600"$(( \
    90 + SLURM_ARRAY_TASK_ID + i * 2))"_2_sequence.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-2600"$(( \
    90 + SLURM_ARRAY_TASK_ID + i * 2))"_2.fastq
done

# Concatenating all relevant microcosm reads into two giant fastq files.
echo "Concatenating all reads"
cat /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-2600*_1.fastq \
    > /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_ALL_READS_1.fastq
cat /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-2600*_2.fastq \
    > /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_ALL_READS_2.fastq

# Removing the symbolic links to avoid overloading.
echo "Removing symbolic links"
rm /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-2600*.fastq

# Deleting unzipped fastq files.
echo "Deleting fastq files"
for ((i=0; i<8; i++)); do
    rm -r \
    /home/akshitg/corderolab/akshitg/D19-2600"$((90 + SLURM_ARRAY_TASK_ID + i * 2))"-4151H/
done

# Running the metaspades assembly on the microcosm reads.
echo "Running metaspades assembly"
metawrap assembly \
    -1 \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_ALL_READS_1.fastq \
    -2 \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_ALL_READS_2.fastq \
    -m 95 \
    -t 20 \
    --metaspades \
    -o /home/akshitg/corderolab/akshitg/Assemblies/M0"${SLURM_ARRAY_TASK_ID}"-assembly

# Removing the combined reads.
rm /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_ALL_READS*.fastq

# Creating a directory called done in the assembly folder.
mkdir /home/akshitg/corderolab/akshitg/Assemblies/M0"${SLURM_ARRAY_TASK_ID}"-assembly/done
