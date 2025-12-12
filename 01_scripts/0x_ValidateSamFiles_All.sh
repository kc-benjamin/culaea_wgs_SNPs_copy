#!/bin/bash

#SBATCH --partition=batch  
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=0-24:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out
#SBATCH --array=0-96

ml Java/17.0.6 GATK/4.6.0.0-GCCcore-13.2.0-Java-17

#Global variables
BAM_FOLDERtest="06_bam_files/test2"
LOG_FOLDER="98_log_files"

#Array number
samp_num=$(($SLURM_ARRAY_TASK_ID +1))
# Fetch filename from the array
sample_name=$(cut -f1 02_info_files/SRR_Acc_List_ML.txt | sed -n "${samp_num}p")

#File validation
gatk ValidateSamFile \
    -I $BAM_FOLDERtest/${sample_name}.dedup.bam \
    -O ${LOG_FOLDER}/test2Val_${SLURM_ARRAY_TASK_ID}.log \
    --MODE SUMMARY