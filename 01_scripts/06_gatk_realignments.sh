#!/bin/bash

#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=10G
#SBATCH --time=00-24:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out
#SBATCH --error=98_log_files/%x_%j.err
#SBATCH --array=0-96

#PREFIX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" 02_info_files/SRR_Acc_List_ML.txt)
# Copy script to log folder
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98_log_files"
cp "$SCRIPT" "$LOG_FOLDER"/"$TIMESTAMP"_"$NAME"

# Load needed modules
module load SAMtools/1.18-GCC-12.3.0

# Global variables
BAM="06_bam_files/test2"
GENOMEFOLDER="03_genome"
GENOME=$(ls -1 $GENOMEFOLDER/brook_genome_hap1_v1.fa | xargs -n 1 basename) #changed to brook genome
GENOME_FULL="$GENOMEFOLDER/$GENOME"

# Build Bam Index
echo " >>> Realigning...
"

#Pass the sample number from the sbatch command
samp_num=$((SLURM_ARRAY_TASK_ID +1))

# Fetch filename from the array
sample_name=$(cut -f1 02_info_files/SRR_Acc_List_ML.txt | sed -n "${samp_num}p")
file=${sample_name}.dedup.bam ##ensure this file exists

echo "
     >>> Realigning TARGET for $file <<<
     "

# Index the bam file First #what is this for?
samtools index $BAM/$file

# Now load modules
module purge
#module load nixpkgs/16.09 #ensure that this is right#
ml Java/1.8.0_241 GATK/3.8-1-Java-1.8.0_241

# Realign
java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R $GENOMEFOLDER/$GENOME \
    -I $BAM/$file \
    -o $BAM/${sample_name}.intervals

echo "
     >>> Realigning INDELs for $file <<<
     "
java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R $GENOMEFOLDER/$GENOME \
    -I $BAM/$file \
    -targetIntervals $BAM/${sample_name}.intervals \
    --consensusDeterminationModel USE_READS  \
    -o $BAM/${sample_name}.realigned.bam

echo ">>> Cleaning a bit...
"
echo "
DONE! Check your files"
