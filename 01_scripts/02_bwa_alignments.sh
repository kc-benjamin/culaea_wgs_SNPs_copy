#!/bin/bash

#SBATCH --partition=batch  
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=0-24:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out

# Load needed modules
module load bwa samtools

# Global variables
GENOMEFOLDER="03_genome"
GENOME=$(ls -1 $GENOMEFOLDER/*fa | xargs -n 1 basename)
INDGENOME=${GENOME}.fai
RAWDATAFOLDER="05_trimmed_data"
ALIGNEDFOLDER="06_bam_files"
LOG_FOLDER="98_log_files"

# Test if user specified a number of CPUs
if [[ -z "$NCPU" ]]
then
    NCPU=4
fi

#Pass the sample number from the sbatch command
samp_num=$1

# Pull sample name from the sample info
name=$(cut -f1 /home/kcb95328/culaea_wgs_SNPs_copy/02_info_files/datatable.txt | sed -n "${samp_num}p")

# Name of uncompressed file
file1=${name}.R1.trimmed.fastq.gz
file2=${name}.R2.trimmed.fastq.gz
echo ">>> Aligning file $file1 $file2 <<<"

# Set read group header line info
RG="@RG\tID:${name}\tSM:${name}\tPL:Illumina"

# Align reads
bwa mem -t $NCPU -R $RG $GENOMEFOLDER/$GENOME $RAWDATAFOLDER/$file1 $RAWDATAFOLDER/$file2 |
samtools view -Sb -q 10 - > $ALIGNEDFOLDER/${name%}.bam

# Sort
samtools sort --threads $NCPU $ALIGNEDFOLDER/${name%.R1.trimmed.fastq.gz}.bam \
    > $ALIGNEDFOLDER/${name%.R1.trimmed.fastq.gz}.sorted.bam

# Index
samtools index $ALIGNEDFOLDER/${name%.R1.trimmed.fastq.gz}.sorted.bam

&> $LOG_FOLDER/02_mapping_${name}.log
