#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=10G
#SBATCH --time=00-24:00:00
#SBATCH --account=def-jonmee
#SBATCH --mail-user=jmee@mtroyal.ca
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out

# Copy script to log folder
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98_log_files"
cp "$SCRIPT" "$LOG_FOLDER"/"$TIMESTAMP"_"$NAME"

# Load needed modules
module load StdEnv/2020 samtools/1.12

# Global variables
BAM="06_bam_files"
GENOMEFOLDER="03_genome"
GENOME=$(ls -1 $GENOMEFOLDER/*fa | xargs -n 1 basename)

# Build Bam Index
echo " >>> Realigning...
"

#Pass the sample number from the sbatch command
samp_num=$1

# Fetch filename from the array
sample_name=$(cut -f1 02_info_files/datatable.txt | sed -n "${samp_num}p")
file=${sample_name}.dedup.bam

echo "
     >>> Realigning TARGET for $file <<<
     "

# Index the bam file First
samtools index $BAM/$file

# Now load modules
module purge
module load nixpkgs/16.09
module load java gatk/3.8

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
