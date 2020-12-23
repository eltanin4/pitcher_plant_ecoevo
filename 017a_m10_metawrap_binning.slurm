#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --array=10
#SBATCH -t 12:00:00
#SBATCH -p sched_mit_chisholm,newnodes,sched_mit_chisholm,sched_mit_hill,sched_any
#SBATCH -J bin_m10
#SBATCH --mem=100G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/metawrap_binning_%a.out

# Adding symbolic links to all fastq files.
echo "Creating symbolic links"
for ((i=0; i<8; i++)); do
    ln -s \
    /home/akshitg/corderolab/akshitg/D19-260"$((90 + SLURM_ARRAY_TASK_ID + i * 2))"-4151H/190801Cor_D19-260"$(( \
    90 + SLURM_ARRAY_TASK_ID + i * 2))"_1_sequence.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260"$(( \
    90 + SLURM_ARRAY_TASK_ID + i * 2))"_1.fastq
    ln -s \
    /home/akshitg/corderolab/akshitg/D19-260"$((90 + SLURM_ARRAY_TASK_ID + i * 2))"-4151H/190801Cor_D19-260"$(( \
    90 + SLURM_ARRAY_TASK_ID + i * 2))"_2_sequence.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260"$(( \
    90 + SLURM_ARRAY_TASK_ID + i * 2))"_2.fastq
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
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260100_1.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260102_1.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260104_1.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260108_1.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260110_1.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260112_1.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260114_1.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260100_2.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260102_2.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260104_2.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260108_2.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260110_2.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260112_2.fastq \
    /home/akshitg/corderolab/akshitg/M0"${SLURM_ARRAY_TASK_ID}"_190801Cor_D19-260114_2.fastq

