#!/bin/zsh
#SBATCH -N 1
#SBATCH -c 20
#SBATCH -t 18:00:00
#SBATCH --array=1-10
#SBATCH -p sched_mit_chisholm
#SBATCH -J virsort
#SBATCH --mem=100G
#SBATCH --output=/home/akshitg/lte/metagenomes_lte/slurm_out/virsort_%a.out

echo "Executing VirSorter run . . ."
echo "Creating and removing old traces and temporary directories."
rm -r /home/akshitg/corderolab/akshitg/Binning/M0"$((SLURM_ARRAY_TASK_ID))"-binning/viroutdir
conda init zsh
conda activate virsorter

echo "Concatenating all assembled bins into one file."
cat \
    /home/akshitg/corderolab/akshitg/Binning/M0"$((SLURM_ARRAY_TASK_ID))"-binning/BIN_REFINEMENT/metawrap_50_10_bins/*.fa \
    > /home/akshitg/corderolab/akshitg/Binning/M0"$((SLURM_ARRAY_TASK_ID))"-binning/BIN_REFINEMENT/metawrap_50_10_bins/assembled_bins.fa

echo "Detecting phages in assembled bins."
~/miniconda2/envs/virsorter/bin/wrapper_phage_contigs_sorter_iPlant.pl \
    -f /home/akshitg/corderolab/akshitg/Binning/M0"$((SLURM_ARRAY_TASK_ID))"-binning/BIN_REFINEMENT/metawrap_50_10_bins/assembled_bins.fa \
    --db 1 \
    --wdir /home/akshitg/corderolab/akshitg/Binning/M0"$((SLURM_ARRAY_TASK_ID))"-binning/viroutdir \
    --ncpu 20 \
    --no_c \
    --keep-db \
    --diamond \
    --data-dir /home/akshitg/virsorter_data_dir/virsorter-data

 echo "Virsorter run complete."
