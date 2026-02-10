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
#SBATCH --error=98_log_files/%x_%j.err
#SBATCH --array=1-23

# Load needed modules
module load BCFtools/1.21-GCC-13.3.0
#gcc/9.3.0

# Global variables
INFO="02_info_files"
GENOMEFOLDER="03_genome"
GENOME=$(ls -1 $GENOMEFOLDER/brook_genome_hap1_v1.fa | xargs -n 1 basename)
VCF="07_raw_VCFs"
BAM="02_info_files/ML_bamfiles.txt"
SAMPS="02_info_files/SRR_Acc_List_ML.txt" #why does this not split it by file?

#Pass the chromosome number from the sbatch command
chrom_num=$1

# Fetch chromosome from the array
CHROM=$(sed -n "${chrom_num}p" 03_genome/brook_genome_hap1_v1_chromosomes2.txt)
SCAFFOLD=$(echo "$CHROM" | grep -oP 'scaffold\d+')


bcftools mpileup -Ou --fasta-ref $GENOMEFOLDER/$GENOME --bam-list "$BAM" -q 5 -r $CHROM -I -a FMT/AD | \
	bcftools call -S "$SAMPS" -G - -f GQ -mv -Ov > "$VCF/${SCAFFOLD}.vcf"
