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
#SBATCH --array=1-119

#PREFIX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" 02_info_files/SRR_Acc_List_ML.txt)
# Load needed modules
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh 
conda activate /home/kcb95328/conda/envs/culaea_pkgs

# Global variables
POPULATION="/home/kcb95328/Info-Astotin" #change this line as needed

GENOMEFOLDER="/home/kcb95328/genomes"
GENOME=$(ls -1 $GENOMEFOLDER/brook_genome_hap1_v1_amhy_masked.fa | xargs -n 1 basename)
GENOME_FULL="$GENOMEFOLDER/$GENOME"

RAWDATAFOLDER="/scratch/kcb95328/AstotinLakeBrooks/05_trimmed_data" #change as needed to the main population directory
ALIGNEDFOLDER="01_aligned_bams" #change as needed; should be from the submit directory
LOG_FOLDER="98_log_files"
#echo "$GENOME and $INDGENOME found in $GENOMEFOLDER"

# Test if user specified a number of CPUs
if [[ -z "$NCPU" ]]
then
    NCPU=4
fi

#Pass the sample number from the sbatch command
samp_num=$SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_TASK_ID='$SLURM_ARRAY_TASK_ID'"
echo "samp_num='$samp_num'"
#echo "PREFIX='$PREFIX'"

# Pull sample name from the sample info
name=$(cut -f1 $POPULATION/SRR_Acc_List_AL.txt | sed -n "${samp_num}p") #change as needed
#name=SRR19221338

# Name of uncompressed file
file1=${name}.R1.trimmed.fastq.gz
file2=${name}.R2.trimmed.fastq.gz
#file1=${name}.fastq
echo ">>> Aligning file $file1 and $file2 <<<"

# Set read group header line info
#RG=$'@RG\tID:'"${name}"$'\tSM:'"${name}"$'\tPL:Illumina'
RG="@RG\tID:${name}\tSM:${name}\tPL:Illumina"
echo $RG

# Align reads
#bwa index $GENOME_FULL bwa-generated-index
bwa mem -t $NCPU -R $RG $GENOME_FULL $RAWDATAFOLDER/$file1 $RAWDATAFOLDER/$file2 |
    samtools view -b -o "$ALIGNEDFOLDER/${name}.bam" #previously used -q 10 to generate first sets of files (or -q 20 for the Increased Filtering files)

# Sort
samtools sort -@ $NCPU $ALIGNEDFOLDER/${name}.bam \
    -o $ALIGNEDFOLDER/${name}.trimmed.fastq.gz.sorted.bam

# Index
samtools index $ALIGNEDFOLDER/${name}.trimmed.fastq.gz.sorted.bam
    &> $LOG_FOLDER/02_mapping_${name}.log

conda deactivate