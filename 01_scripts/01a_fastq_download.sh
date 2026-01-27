#!/bin/bash

#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=0-12:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%a_%x.out
#SBATCH --error=98_log_files/%a_%x.err
#SBATCH --array=1-97

#Go to correct directory
cd $SLURM_SUBMIT_DIR
#load package
ml SRA-Toolkit/3.0.3-gompi-2022a
#get accession number
name=$(sed -n "$(($SLURM_ARRAY_TASK_ID))p" /scratch/kcb95328/Mee-Culaea-WGS/02_info_files/SRR_Acc_List_ML.txt)
#fastq conversion
fasterq-dump $name --split-spot