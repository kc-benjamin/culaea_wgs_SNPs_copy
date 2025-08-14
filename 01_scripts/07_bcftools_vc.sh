#!/bin/bash

#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=1024M
#SBATCH --time=0-12:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out

# Load needed modules
module gcc/9.3.0 bcftools/1.11

# Global variables
INFO="02_info_files"
GENOMEFOLDER="03_genome"
GENOME=$(ls -1 $GENOMEFOLDER/*fa | xargs -n 1 basename)
VCF="./07_raw_VCFs"
BAM="./02_info_files/bamfiles.txt"
SAMPS="./02_info_files/sample_file.txt"

#Pass the chromosome number from the sbatch command
chrom_num=$1

# Fetch chromosome from the array
CHROM=$(sed -n "${chrom_num}p" ./02_info_files/chromosomes.txt)

bcftools mpileup -Ou -f $GENOMEFOLDER/$GENOME --bam-list $BAM -q 5 -r $CHROM -I -a FMT/AD | bcftools call -S $SAMPS -G - -f GQ -mv -Ov > $VCF/${CHROM}\.vcf
