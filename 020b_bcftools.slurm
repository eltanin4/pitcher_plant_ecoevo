#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --array=108
#SBATCH -t 12:00:00
#SBATCH -p newnodes,sched_mit_chisholm,sched_mit_hill,sched_any
#SBATCH -J bcfmg
#SBATCH --mem=50G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/bcftools_mg_%a.out

bcftools mpileup \
    --threads 10 \
    -a FORMAT/AD \
    -Q 15 \
    -L 1000 \
    -d 100000 \
    -m 2 \
    -f /home/akshitg/lte/sags_lte/combined.fasta \
    /home/akshitg/lte/metagenomes_lte/alignments/D19-260"$((SLURM_ARRAY_TASK_ID))"_minimap2_axsr_results.bam | \
    bcftools call \
    --threads 10 \
    --ploidy 1 \
    -Ov -mv \
    -o /home/akshitg/lte/metagenomes_lte/freebayes_out/mg_"$((SLURM_ARRAY_TASK_ID))"_combined_bcftools.vcf
