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

module load SAMtools/1.21-GCC-13.3.0

GENOMEFOLDER="03_genome"
GENOME=$(ls -1 $GENOMEFOLDER/brook_genome_hap1_v1.fa | xargs -n 1 basename) #changed to brook genome
GENOME_FULL="$GENOMEFOLDER/$GENOME"
##INDGENOME="${GENOME}.fai"
#RAWDATAFOLDER="05_trimmed_data"
#ALIGNEDFOLDER="06_bam_files"
#ALIGNED_test="06_bam_files/test3"
LOG_FOLDER="98_log_files"
#echo "$GENOME and $INDGENOME found in $GENOMEFOLDER"
name=$(cut -f1 02_info_files/SRR_Acc_List_ML.txt | sed -n "${samp_num}p")

samtools view -q 10 -o "06_bam_files/${name}.sam" 06_bam_files/${name}.bam