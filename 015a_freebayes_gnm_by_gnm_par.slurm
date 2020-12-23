#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --array=35-99
#SBATCH -t 3:00:00
#SBATCH -p newnodes,sched_mit_chisholm,sched_mit_hill,sched_any
#SBATCH -J fbg1
#SBATCH --mem=30G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/fb_par_gnm_1_mg_%a.out

/home/akshitg/.linuxbrew/bin/freebayes-parallel /home/akshitg/lte/metagenomes_lte/freebayes_out/gnm_1_regions.txt 10 \
        -f /home/akshitg/lte/sags_lte/combined.fasta -p 1 \
        --haplotype-length 0 --min-alternate-count 1 \
        --min-alternate-fraction 0 --pooled-continuous \
        --theta 0.001 --ploidy 1 -J -K -X -u -n 0 --haplotype-length 3 --min-repeat-size 5 \
        --min-repeat-entropy 0 -m 20 -q 30 -R 0 -Y 0 -e 1000 -F 0.2 -C 2 -G 1 \
        --min-coverage 0 --min-alternate-qsum 0 \
        /home/akshitg/lte/metagenomes_lte/D19-2600"$((SLURM_ARRAY_TASK_ID))"_minimap2_axsr_results.bam \
        > /home/akshitg/lte/metagenomes_lte/freebayes_out/mg_"$((SLURM_ARRAY_TASK_ID))"_gnm_1_only_snps_par.vcf

