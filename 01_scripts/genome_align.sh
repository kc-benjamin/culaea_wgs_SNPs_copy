#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name="genome_align"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00-12:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out
#SBATCH --error=98_log_files/%x_%j.err

module load MUMmer/4.0.1-GCCcore-13.3.0

nucmer -t 4 -c 100 -p attempt2 stickleback_v5_assembly.fa brook_genome_hap1_v1.fa

show-coords attempt2.delta > coords.txt

mummerplot --medium -t png attempt2.delta