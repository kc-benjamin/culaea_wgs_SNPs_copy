#!/bin/bash

#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=0-06:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out

# Load modules
module java/13.0.2 picard/2.26.3

# Global variables
PICARD=$EBROOTPICARD/picard.jar
MARKDUPS="MarkDuplicates"
ALIGNEDFOLDER="06_bam_files"
METRICSFOLDER="99_metrics"

# Copy script to log folder
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98_log_files"

export JAVA_TOOL_OPTIONS="-Xms2g -Xmx50g "
export _JAVA_OPTIONS="-Xms2g -Xmx50g "

#Pass the sample number from the sbatch command
samp_num=$1

# Fetch filename from the array
sample_name=$(cut -f1 02_info_files/datatable.txt | sed -n "${samp_num}p")
file=${sample_name}.sorted.bam

echo "DEduplicatING sample $file"

java -jar $PICARD $MARKDUPS \
    INPUT=$ALIGNEDFOLDER/$file \
    OUTPUT=$ALIGNEDFOLDER/${sample_name}.dedup.bam \
    METRICS_FILE=$METRICSFOLDER/${sample_name}_DUP_metrics.txt \
    VALIDATION_STRINGENCY=SILENT \
    REMOVE_DUPLICATES=true

echo " >>> Cleaning a bit...
"
echo "DONE! Go check your files."
