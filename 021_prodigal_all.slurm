#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 10
#SBATCH -t 12:00:00
#SBATCH -p newnodes,sched_mit_chisholm,sched_mit_hill,sched_any
#SBATCH -J prod
#SBATCH --mem=50G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/prodigal_all.out

for ((i=1; i<=9; i++));
   do
        prodigal \
            -i /home/akshitg/main/lte/sags_lte/D19-26000${i}/contigs_gnm.fasta \
            -o /home/akshitg/main/lte/sags_lte/D19-26000${i}/gene_list.faa \
            -d /home/akshitg/main/lte/sags_lte/D19-26000${i}/genes.faa \
            -a /home/akshitg/main/lte/sags_lte/D19-26000${i}/proteins.faa
   done
for ((i=10; i<=24; i++));
    do
        prodigal \
            -i /home/akshitg/main/lte/sags_lte/D19-2600${i}/contigs_gnm.fasta \
            -o /home/akshitg/main/lte/sags_lte/D19-2600${i}/gene_list.faa \
            -d /home/akshitg/main/lte/sags_lte/D19-2600${i}/genes.faa \
            -a /home/akshitg/main/lte/sags_lte/D19-2600${i}/proteins.faa
    done
for ((i=26; i<=34; i++));
    do
        prodigal \
            -i /home/akshitg/main/lte/sags_lte/D19-2600${i}/contigs_gnm.fasta \
            -o /home/akshitg/main/lte/sags_lte/D19-2600${i}/gene_list.faa \
            -d /home/akshitg/main/lte/sags_lte/D19-2600${i}/genes.faa \
            -a /home/akshitg/main/lte/sags_lte/D19-2600${i}/proteins.faa
    done
