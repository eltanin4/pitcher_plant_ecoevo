#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH --array=10-34
#SBATCH -p sched_mit_chisholm
#SBATCH -J fbnpm01
#SBATCH --mem=30G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/fbnpm01_%j_%a.out

# mkdir /home/akshitg/lte/metagenomes_lte/freebayes_tmp
export TMPDIR=/home/akshitg/lte/metagenomes_lte/freebayes_tmp
# mkdir /home/akshitg/lte/metagenomes_lte/freebayes_out

echo Running freebayes non-parallel for M01, array "$((SLURM_ARRAY_TASK_ID))"
freebayes -f /home/akshitg/main/lte/sags_lte/D19-2600"$((SLURM_ARRAY_TASK_ID))"/contigs_gnm.fasta -p 1 \
        --haplotype-length 0 --min-alternate-count 1 \
        --min-alternate-fraction 0 --pooled-continuous \
        --report-monomorphic \
           /home/akshitg/lte/metagenomes_lte/D19-260035_minimap2_axsr_results.bam \
        > /home/akshitg/lte/metagenomes_lte/freebayes_out/m01_35_gnm_"$((SLURM_ARRAY_TASK_ID))"_freebayes_naive.vcf

