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

module load SAMtools/1.21-GCC-13.3.0

samp_num=$(($SLURM_ARRAY_TASK_ID +1))
name=$(cut -f1 02_info_files/SRR_Acc_List_ML.txt | sed -n "${samp_num}p")

samtools sort -@ $NCPU 04_trimmed_data/${name}.fastq  #\
    #-o 06_bam_files/${name}.R1.trimmed.fastq.gz.bam