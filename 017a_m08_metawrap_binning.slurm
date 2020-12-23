#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --array=8
#SBATCH -t 12:00:00
#SBATCH -p sched_mit_chisholm,newnodes,sched_mit_chisholm,sched_mit_hill,sched_any
#SBATCH -J bin_m08
#SBATCH --mem=100G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/metawrap_binning_%a.out

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

# Making the binning and initial directories
rm -r /home/akshitg/corderolab/akshitg/Binning/M0"${SLURM_ARRAY_TASK_ID}"-binning/
mkdir /home/akshitg/corderolab/akshitg/Binning/M0"${SLURM_ARRAY_TASK_ID}"-binning
mkdir /home/akshitg/corderolab/akshitg/Binning/M0"${SLURM_ARRAY_TASK_ID}"-binning/INITIAL_BINNING
mkdir /home/akshitg/corderolab/akshitg/Binning/M0"${SLURM_ARRAY_TASK_ID}"-binning/INITIAL_BINNING/work_files

# Running the binner.
metaWRAP binning \
    -o /home/akshitg/corderolab/akshitg/Binning/M0"${SLURM_ARRAY_TASK_ID}"-binning/INITIAL_BINNING \
    -t 20 \
    -m 95 \
    -a /home/akshitg/corderolab/akshitg/Assemblies/M0"${SLURM_ARRAY_TASK_ID}"-assembly/out/final.contigs.fa \
    --metabat2 --maxbin2 --concoct \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260042_1.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260050_1.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260058_1.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260066_1.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260074_1.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260082_1.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260090_1.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260098_1.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260042_2.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260050_2.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260058_2.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260066_2.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260074_2.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260082_2.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260090_2.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260098_2.fastq

