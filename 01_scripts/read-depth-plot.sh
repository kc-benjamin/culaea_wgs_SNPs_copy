#!/bin/bash

#SBATCH --partition=highmem_p
#SBATCH --job-name="read_depth_R"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G
#SBATCH --time=0-24:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch/kcb95328/ShundaLakeBrooks/98_log_files/read_depth_R_%j.out
#SBATCH --error=/scratch/kcb95328/ShundaLakeBrooks/98_log_files/read_depth_R_%j.err

#go to submit directory
cd $SLURM_SUBMIT_DIR

#activate R script
ml R/4.5.1-gfbf-2025a

R --no-save < /home/kcb95328/culaea_wgs_SNPs_copy/01_scripts/read-depth-plot.r
