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
#SBATCH --array=1-1

#PREFIX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" 02_info_files/SRR_Acc_List_ML.txt)
# Load needed modules
module load BWA/0.7.18-GCCcore-13.3.0 SAMtools/1.21-GCC-13.3.0

# Global variables
GENOMEFOLDER="/scratch/kcb95328/Mee-Culaea-WGS/03_genome"
GENOME=$(ls -1 $GENOMEFOLDER/brook_genome_hap1_v1.fa | xargs -n 1 basename) #changed to ninespine genome
GENOME_FULL="$GENOMEFOLDER/$GENOME"
##INDGENOME="${GENOME}.fai"
RAWDATAFOLDER="05_trimmed_data"
ALIGNEDFOLDER="06_bam_files"
#ALIGNED_test="06_bam_files/test3"
LOG_FOLDER="98_log_files"
echo "$GENOME and $INDGENOME found in $GENOMEFOLDER"

# Test if user specified a number of CPUs
if [[ -z "$NCPU" ]]
then
    NCPU=4
fi

#Pass the sample number from the sbatch command
samp_num=$SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_TASK_ID='$SLURM_ARRAY_TASK_ID'"
echo "samp_num='$samp_num'"
#echo "PREFIX='$PREFIX'"

# Pull sample name from the sample info
name=$(cut -f1 02_info_files/SRR_Acc_List_ML.txt | sed -n "${samp_num}p")

# Name of uncompressed file
file1=${name}.R1.trimmed.fastq.gz
file2=${name}.R2.trimmed.fastq.gz
#file1=${name}.fastq
echo ">>> Aligning file $file1 and $file2 <<<"

# Set read group header line info
#RG=$'@RG\tID:'"${name}"$'\tSM:'"${name}"$'\tPL:Illumina'
RG="@RG\tID:${name}\tSM:${name}\tPL:Illumina"
echo $RG

# Align reads
#bwa index $GENOME_FULL bwa-generated-index
bwa mem -t $NCPU -R $RG $GENOME_FULL $RAWDATAFOLDER/$file1 $RAWDATAFOLDER/$file2 |
    samtools view -b -q 10 -o "$ALIGNEDFOLDER/${name}.bam" 

# Sort
samtools sort -@ $NCPU $ALIGNEDFOLDER/${name}.bam \
    -o $ALIGNEDFOLDER/${name}.trimmed.fastq.gz.sorted.bam

# Index
#samtools index $ALIGNEDFOLDER/${name}.trimmed.fastq.gz.sorted.bam
    #&> $LOG_FOLDER/02_mapping_${name}.log
