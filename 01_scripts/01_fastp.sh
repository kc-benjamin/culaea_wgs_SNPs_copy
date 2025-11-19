#!/bin/bash

#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=10G
#SBATCH --time=0-00:30:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%a_%x.out
#SBATCH --error=98_log_files/%a_%x.err
#SBATCH --array=0-96

PREFIX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" 02_info_files/SRR_Acc_List_ML.txt)
#currently for Muir Lake
# Load up fastp
module load fastp/0.23.4-GCC-13.2.0

# Variables
INDIR="/scratch/kcb95328/Mee-Culaea-WGS/04_raw_data"
OUTDIR="/scratch/kcb95328/Mee-Culaea-WGS/05_trimmed_data"
LOG="/scratch/kcb95328/Mee-Culaea-WGS/98_log_files"
#mkdir $OUTDIR/01_reports
#need to make the directory 01_reports in 05_trimmed_data

#Pass the sample number from the sbatch command
#samp_num=$PREFIX

#Pull sample name from the sample info
sample_name=$(cut -f1 02_info_files/SRR_Acc_List_ML.txt | sed -n "${$PREFIX}p")

fastp -w ${SLURM_CPUS_PER_TASK} \
        -i "$INDIR/$(cut -f12 02_info_files/SRR_Acc_List_ML.txt | sed -n "${$PREFIX}p").fastq" \
        -I "$INDIR/$(cut -f13 02_info_files/SRR_Acc_List_ML.txt | sed -n "${$PREFIX}p").fastq" \
        -o $OUTDIR/"$sample_name".R1.trimmed.fastq.gz \
        -O $OUTDIR/"$sample_name".R2.trimmed.fastq.gz \
        -j $OUTDIR/01_reports/"$sample_name".json \
        -h $OUTDIR/01_reports/"$sample_name".html \
        &> "$LOG"/01_fastp_"$sample_name".out
