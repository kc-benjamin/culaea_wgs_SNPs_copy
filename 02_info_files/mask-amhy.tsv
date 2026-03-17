#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name="mask-amhy"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00-12:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out
#SBATCH --error=98_log_files/%x_%j.err

ml BEDTools/2.31.1-GCC-13.3.0 SAMtools/1.21-GCC-13.3.0

### mask.bed file < made with nano
### PGA_scaffold14__88_contigs__length_21401847 15,445,433  15,450,428

### mask Y par in reference fasta
GENOME="03_genome"
bedtools maskfasta \
    -bed $GENOME/mask-amhy.bed \
    -fi $GENOME/brook_genome_hap1_v1.fa \
    -fo $GENOME/brook_genome_hap1_v1_no_amhy.fa

### index masked fasta
samtools faidx $GENOME/brook_genome_hap1_v1_no_amhy.fa