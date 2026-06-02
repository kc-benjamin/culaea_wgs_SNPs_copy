#!/bin/bash
#SBATCH --partition=batch
#SBATCH --partition=batch
#SBATCH --job-name="star-align"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00
#SBATCH --mem=1G
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=star-align_%j.out
#SBATCH --error=star-align_%j.err
#SBATCH --array=1-40

ml STAR/2.7.11a-GCC-12.3.0

#variables
INDIR="01_raw_reads"
OUTDIR="03_star_output"
GENOMEDIR="02_genome"

#Pass the sample number from the sbatch command
samp_num=$SLURM_ARRAY_TASK_ID

# Pull sample name from the sample info
name=$(cut -f1 raw_data_list3.txt | sed -n "${samp_num}p")

STAR --runThreadN 4 \
     --genomeDir $GENOMEDIR \
     --readFilesIn $INDIR/${name}_R1_001.fastq.gz $INDIR/${name}_R2_001.fastq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTDIR/${name} \
     --outSAMtype BAM SortedByCoordinate

