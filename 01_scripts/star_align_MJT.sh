#!/bin/bash
#SBATCH --partition=batch          # use other queues for >7 days, more memory, gpu
#SBATCH --job-name=star_align
#SBATCH --ntasks=1
#SBATCH --time=24:00:00         # D-HH:MM:SS format
#SBATCH --mem=100G                # memory per node
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=99_log_files/star-align_%j.out
#SBATCH --error=99_log_files/star-align_%j.err
#SBATCH --array=1-40                # set range

set -euo pipefail
date | tee /dev/stderr
cd $SLURM_SUBMIT_DIR

ml STAR/2.7.11b-GCC-13.3.0

SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" raw_data_list3.txt)
INDIR="/scratch/kcb95328/brook_RNA_Mar2025/01_raw_reads"
R1="${INDIR}/${SAMPLE}_R1_001.fastq.gz"
R2="${INDIR}/${SAMPLE}_R2_001.fastq.gz"
OUTDIR="/scratch/kcb95328/brook_RNA_Mar2025/04_counts"
INDEX="/scratch/kcb95328/brook_RNA_Mar2025/02_genome/star_index"


# readFilesCommand only needed if reads are zipped (.gz)
STAR --runMode alignReads \
    --runThreadN 8 \
    --genomeDir "${INDEX}" \
    --readFilesIn "${R1}" "${R2}" \
    --readFilesCommand zcat \
    --quantMode GeneCounts \
    --outFilterMultimapNmax 1 \
    --outFileNamePrefix "${OUTDIR}/${SAMPLE}_" \
    --outSAMtype None




### lines to update
# 8 your email address
# 10 & 11 - log file locations
# 12 - number of jobs = number of samples
# 15 - target directory for log files
# 19 - starting location that contains samples.txt
# 25 - directory of sequence files
# 28 - output directory for results
# 29 - location of reference genome index files
