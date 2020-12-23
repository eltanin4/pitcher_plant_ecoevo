#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH -p sched_mit_chisholm, sched_mit_hill, sched_any, newnodes
#SBATCH -J cp_mm2_100
#SBATCH --mem=100G

for ((i=100; i<=114; i++));
    do
        cat "~/lte/metagenomes_lte/D19-260"$i"-4151H/190801Cor_D19-260"$i"_1_sequence.fastq.gz" "~/lte/metagenomes_lte/D19-260"$i"-4151H/190801Cor_D19-260"$i"_1_sequence.fastq.gz" >> "metagenomes_lte/D19-260"$i"-4151H_combined_mg_reads.fastq.gz"
        minimap2 ~/lte/sags_lte/combined.fasta "~/lte/metagenomes_lte/D19-260"$i"-4151H_combined_mg_reads.fastq.gz" > "~/lte/metagenomes_lte/D19-260"$i"_minimap2_default_results.paf"
    done



