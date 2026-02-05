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
#SBATCH --array=1-97

module load FastQC/0.12.1-Java-11

fastqc -o /scratch/kcb95328/Mee-Culaea-WGS/05_trimmed_data/05_fastQC-new-split-files \
    -t ${SLURM_CPUS_PER_TASK} \
    05_trimmed_data/$(sed -n "${SLURM_ARRAY_TASK_ID}p" 02_info_files/SRR_Acc_List_ML.txt)R1.trimmed.fastq.gz \
    05_trimmed_data/$(sed -n "${SLURM_ARRAY_TASK_ID}p" 02_info_files/SRR_Acc_List_ML.txt)R2.trimmed.fastq.gz