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
#SBATCH --array=1-1

# Load needed modules
module load BCFtools/1.21-GCC-13.3.0
#gcc/9.3.0

# Global variables
INFO="02_info_files"
GENOMEFOLDER="03_genome"
GENOME=$(ls -1 $GENOMEFOLDER/brook_genome_hap1_v1.fa | xargs -n 1 basename)
VCF="07_raw_VCFs"
BAM="02_info_files/ML_bamfiles_full.txt"
SAMPS="02_info_files/SRR_Acc_List_ML.txt" #why does this not split it by file?

#Pass the chromosome number from the sbatch command
chrom_num=$1

# Fetch chromosome from the array
CHROM=$(sed -n "${chrom_num}p" 03_genome/brook_genome_hap1_v1_chromosomes2.txt)
SCAFFOLD=$(echo "$CHROM" | grep -oP 'scaffold\d+')

# Filter the BAM list to only files that actually exist (github)
BAM_FILTER=$(mktemp)
while IFS= read -r bamfile; do
	[ -z "$bamfile" ] && continue
	if [ -f "$bamfile" ]; then
		echo "$bamfile" >> "$BAM_FILTER"
	fi
done < "$BAM"

if [ ! -s "$BAM_FILTER" ]; then
	echo "No existing BAM files found from $BAM; exiting."
	rm -f "$BAM_FILTER"
	exit 0
fi

# Create a filtered sample list containing only samples present in the filtered BAM list
SAMP_FILTER=$(mktemp)
while IFS= read -r sample; do
	[ -z "$sample" ] && continue
	if grep -Fq "$sample" "$BAM_FILTER"; then
		echo "$sample" >> "$SAMP_FILTER"
	fi
done < "$SAMPS"

if [ ! -s "$SAMP_FILTER" ]; then
	echo "No matching samples found in $BAM_FILTER; exiting without running bcftools."
	rm -f "$BAM_FILTER" "$SAMP_FILTER"
	exit 0
fi

bcftools mpileup -Ou -f $GENOMEFOLDER/$GENOME --bam-list "$BAM_FILTER" -q 5 -r $CHROM -I -a FMT/AD | \
	bcftools call -S "$SAMP_FILTER" -G - -f GQ -mv -Ov > "$VCF/${SCAFFOLD}.vcf"

# cleanup
rm -f "$SAMP_FILTER" "$BAM_FILTER"
