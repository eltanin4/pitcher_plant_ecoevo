#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --array=4-5
#SBATCH -t 18:00:00
#SBATCH -p sched_mit_chisholm
#SBATCH -J megahit
#SBATCH --mem=150G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/megahit_%a.out

# SLURM_ARRAY_TASK_ID corresponds to microcosm number.
# NOTE: M09 and M10 have to be dealt with differently 
# because of their numbering schemes.
# SLURM_ARRAY_TASK_ID=1

# Adding symbolic links to all fastq files.
echo "Creating symbolic links"
for ((i=0; i<8; i++)); do
    ln -s \
    /home/akshitg/corderolab/akshitg/D19-2600"$((34 + SLURM_ARRAY_TASK_ID + i * 8))"-4151H/190801Cor_D19-2600"$(( \
    34 + SLURM_ARRAY_TASK_ID + i * 8))"_1_sequence.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-2600"$(( \
    34 + SLURM_ARRAY_TASK_ID + i * 8))"_1.fastq
    ln -s \
    /home/akshitg/corderolab/akshitg/D19-2600"$((34 + SLURM_ARRAY_TASK_ID + i * 8))"-4151H/190801Cor_D19-2600"$(( \
    34 + SLURM_ARRAY_TASK_ID + i * 8))"_2_sequence.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-2600"$(( \
    34 + SLURM_ARRAY_TASK_ID + i * 8))"_2.fastq
done

# Concatenating all relevant microcosm reads into two giant fastq files.
echo "Concatenating all reads"
cat /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-2600*_1.fastq \
    > /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_ALL_READS_1.fastq
cat /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-2600*_2.fastq \
    > /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_ALL_READS_2.fastq

# Running the megahit assembly on the microcosm reads.
echo "Running megahit assembly"
rm -r /home/akshitg/corderolab/akshitg/Assemblies/M0"${SLURM_ARRAY_TASK_ID}"-assembly/
mkdir /home/akshitg/corderolab/akshitg/Assemblies/M0"${SLURM_ARRAY_TASK_ID}"-assembly

/home/akshitg/megahit/build/megahit \
    -1 \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_ALL_READS_1.fastq \
    -2 \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_ALL_READS_2.fastq \
    -t 20 \
    --min-count 3 \
    --k-min 27 \
    --k-max 127 \
    --k-step 10 \
    -o /home/akshitg/corderolab/akshitg/Assemblies/M0"${SLURM_ARRAY_TASK_ID}"-assembly/out

# Removing the symbolic links to avoid overloading.
echo "Removing symbolic links"
rm /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_ALL_READS_*.fastq

# Creating a directory called done in the assembly folder.
mkdir /home/akshitg/corderolab/akshitg/Assemblies/M0"${SLURM_ARRAY_TASK_ID}"-assembly/done
