#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name="mask-amhy"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00-12:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out
#SBATCH --error=98_log_files/%x_%j.err

module load BEDTools/2.31.1-GCC-13.3.0 SAMtools/1.21-GCC-13.3.0

#mask.bed file contents:
#PGA_scaffold14__88_contigs__length_21401847   15445433       15450428

### mask Y par in reference fasta
bedtools maskfasta \
    -bed 03_genome/mask.bed \
    -fi 03_genome/brook_genome_hap1_v1.fa \
    -fo 03_genome/brook_genome_hap1_v1_amhy_masked.fa

### index masked fasta
samtools faidx 03_genome/brook_genome_hap1_v1_amhy_masked.fa