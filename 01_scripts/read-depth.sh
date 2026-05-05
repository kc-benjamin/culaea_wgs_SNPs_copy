#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name="read-depth-Shunda"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00-12:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out
#SBATCH --error=98_log_files/%x_%j.err

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh 
conda activate /home/kcb95328/conda/envs/culaea_pkgs

#location-specific to exon 7
samtools depth -a -H -J -r "PGA_scaffold14__88_contigs__length_21401847:15,449,686-15,450,142" -f 02_info_files/SL_bamfiles_full.txt -o Shunda-exon7-amhy-depth.txt

conda deactivate