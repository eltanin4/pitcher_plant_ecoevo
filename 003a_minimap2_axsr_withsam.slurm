#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --array=100-114
#SBATCH -t 24:00:00
#SBATCH -p newnodes,sched_mit_chisholm,sched_mit_hill,sched_any
#SBATCH -J 100_mm2
#SBATCH --mem=50G

~/minimap2-2.17_x64-linux/minimap2 \
  -ax sr \
  ~/lte/sags_lte/combined.fasta \
  "$HOME/lte/metagenomes_lte/D19-260"$SLURM_ARRAY_TASK_ID"-4151H/190801Cor_D19-260"$SLURM_ARRAY_TASK_ID"_1_sequence.fastq.gz" \
  "$HOME/lte/metagenomes_lte/D19-260"$SLURM_ARRAY_TASK_ID"-4151H/190801Cor_D19-260"$SLURM_ARRAY_TASK_ID"_2_sequence.fastq.gz" \
  > "D19-260"$SLURM_ARRAY_TASK_ID"_minimap2_axsr_results.sam"

samtools sort "D19-260"$SLURM_ARRAY_TASK_ID"_minimap2_axsr_results.sam" \
  > "D19-260"$SLURM_ARRAY_TASK_ID"_minimap2_axsr_results.bam"

samtools index "D19-260"$SLURM_ARRAY_TASK_ID"_minimap2_axsr_results.bam"
