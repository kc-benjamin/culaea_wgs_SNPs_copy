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

module load BWA/0.7.18-GCCcore-13.3.0

#genome variables
GENOMEFOLDER="/scratch/kcb95328/Mee-Culaea-WGS/03_genome"
GENOME=$(ls -1 $GENOMEFOLDER/brook_genome_hap1_v1.fa | xargs -n 1 basename) #changed to brook genome
GENOME_FULL="$GENOMEFOLDER/$GENOME"
#file locations
RAWDATAFOLDER="05_trimmed_data"
#name=$(cut -f1 02_info_files/SRR_Acc_List_ML.txt | sed -n "${samp_num}p")
# Name of uncompressed file
#file1=${name}.R1.trimmed.fastq.gz
#file2=${name}.R2.trimmed.fastq.gz

bwa mem -M -t $SLURM_CPUS_PER_TASK -R @RG\tID:SRR19221208\tSM:SRR19221208\tPL:Illumina $GENOME_FULL $RAWDATAFOLDER/SRR19221208.R1.trimmed.fastq.gz $RAWDATAFOLDER/SRR19221208.R2.trimmed.fastq.gz