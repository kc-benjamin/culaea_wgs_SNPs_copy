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

module load SRA-Toolkit/3.0.3-gompi-2022a

prefetch GCF_964276395.1 ##note: this is not the UGA assembly, but the NCBI re-analyzed assembly
fastq-dump GCF_964276395.1