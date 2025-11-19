#!/bin/bash

#SBATCH --partition=batch  
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=0-24:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out
#SBATCH --array=0-96

PREFIX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" 02_info_files/SRR_Acc_List_ML.txt)
# Load needed modules
module load BWA/0.7.18-GCCcore-13.3.0 SAMtools/1.18-GCC-12.3.0

# Global variables
GENOMEFOLDER="/scratch/kcb95328/Mee-Culaea-WGS/03_genome"
GENOME=$(ls -1 $GENOMEFOLDER/brook_genome_hap1_v1.fa | xargs -n 1 basename) #changed to brook genome
GENOME_FULL="$GENOMEFOLDER/$GENOME"
INDGENOME="${GENOME_FULL}.fai"
RAWDATAFOLDER="05_trimmed_data"
ALIGNEDFOLDER="06_bam_files"
LOG_FOLDER="98_log_files"
echo "$GENOME and $INDGENOME found in $GENOMEFOLDER"

# Test if user specified a number of CPUs
if [[ -z "$NCPU" ]]
then
    NCPU=4
fi

#Pass the sample number from the sbatch command
samp_num=$((SLURM_ARRAY_TASK_ID + 1))
echo "SLURM_ARRAY_TASK_ID='$SLURM_ARRAY_TASK_ID'"
echo "samp_num='$samp_num'"
echo "PREFIX='$PREFIX'"

# Pull sample name from the sample info
name=$(cut -f1 /scratch/kcb95328/Mee-Culaea-WGS/02_info_files/SRR_Acc_List_ML.txt | sed -n "${samp_num}p")

# Name of uncompressed file
file1=${name}.R1.trimmed.fastq.gz
file2=${name}.R2.trimmed.fastq.gz
echo ">>> Aligning file $file1 $file2 <<<"

# Set read group header line info
RG="@RG\tID:${name}\tSM:${name}\tPL:Illumina"

# Align reads
bwa mem -t $NCPU -R $RG $GENOMEFOLDER/$GENOME $RAWDATAFOLDER/$file1 $RAWDATAFOLDER/$file2 |
samtools view -Sb -q 10 - > $ALIGNEDFOLDER/${name%}.bam
#.........from copilot.......
# check BWA index files; build if missing
MISSING=0
for ext in bwt pac ann amb sa; do
    if [[ ! -f "${GENOME_FULL}.${ext}" ]]; then
     echo "Missing BWA index file: ${GENOME_FULL}.${ext}" >&2
     MISSING=1
   fi
 done
 if [[ $MISSING -eq 1 ]]; then
   echo "Building BWA index for $GENOME_FULL" >&2
   bwa index "$GENOME_FULL" || { echo "bwa index failed" >&2; exit 1; }
 fi
# ensure samtools faidx exists
 if [[ ! -f "${GENOME_FULL}.fai" ]]; then
   echo "Creating samtools faidx for $GENOME_FULL" >&2
   samtools faidx "$GENOME_FULL" || { echo "samtools faidx failed" >&2; exit 1; }
 fi
#......end copilot help.......
# Sort
samtools sort --threads $NCPU $ALIGNEDFOLDER/${name%.R1.trimmed.fastq.gz}.bam \
    > $ALIGNEDFOLDER/${name%.R1.trimmed.fastq.gz}.sorted.bam

# Index
samtools index $ALIGNEDFOLDER/${name%.R1.trimmed.fastq.gz}.sorted.bam

&> $LOG_FOLDER/02_mapping_${name}.log
