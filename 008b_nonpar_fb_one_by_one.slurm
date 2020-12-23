#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 48:00:00
#SBATCH --array=0-7
#SBATCH -p sched_mit_chisholm
#SBATCH -J fbnpm01
#SBATCH --mem=100G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/fbnpm01_%j.out

# mkdir /home/akshitg/lte/metagenomes_lte/freebayes_tmp
export TMPDIR=/home/akshitg/lte/metagenomes_lte/freebayes_tmp
# mkdir /home/akshitg/lte/metagenomes_lte/freebayes_out

echo Running freebayes non-parallel for M01, array "$((SLURM_ARRAY_TASK_ID))"
freebayes -f /home/akshitg/lte/sags_lte/combined.fasta -p 1 \
        --haplotype-length 0 --min-alternate-count 1 \
        --min-alternate-fraction 0 --pooled-continuous \
        --report-monomorphic \
           /home/akshitg/lte/metagenomes_lte/D19-2600"$((35 + SLURM_ARRAY_TASK_ID * 8))"_minimap2_axsr_results.bam \
        > /home/akshitg/lte/metagenomes_lte/freebayes_out/m01_"$((35 + SLURM_ARRAY_TASK_ID * 8))"_freebayes_naive.vcf

#freebayes-parallel <(fasta_generate_regions.py /home/akshitg/lte/sags_lte/combined.fasta.fai 100000) 20 \
#        -f /home/akshitg/lte/sags_lte/combined.fasta -p 1 \
#        --min-alternate-fraction 0 \
#        --report-monomorphic \
#           /home/akshitg/lte/metagenomes_lte/D19-260035_minimap2_axsr_results.bam \
#           /home/akshitg/lte/metagenomes_lte/D19-260043_minimap2_axsr_results.bam \
#           /home/akshitg/lte/metagenomes_lte/D19-260051_minimap2_axsr_results.bam \
#           /home/akshitg/lte/metagenomes_lte/D19-260059_minimap2_axsr_results.bam \
#           /home/akshitg/lte/metagenomes_lte/D19-260067_minimap2_axsr_results.bam \
#           /home/akshitg/lte/metagenomes_lte/D19-260075_minimap2_axsr_results.bam \
#           /home/akshitg/lte/metagenomes_lte/D19-260083_minimap2_axsr_results.bam \
#           /home/akshitg/lte/metagenomes_lte/D19-260091_minimap2_axsr_results.bam \
#        > /home/akshitg/lte/metagenomes_lte/m01_freebayes_parallel.vcf
