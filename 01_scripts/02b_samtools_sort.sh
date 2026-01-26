#!/bin/bash

#SBATCH --partition=batch  
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=6G
#SBATCH --time=0-24:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j_.out
#SBATCH --error=98_log_files/%x_%j_.err

module load SAMtools/1.21-GCC-13.3.0

samtools sort -@ $NCPU 06_bam_files/SRR19221208.sam #\
    #-o 06_bam_files/SRR19221208.trimmed.fastq.gz.sorted.bam