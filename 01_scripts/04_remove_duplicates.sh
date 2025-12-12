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
#SBATCH --error=98_log_files/%x_%j.err
#SBATCH --array=0-96

#PREFIX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" 02_info_files/SRR_Acc_List_ML.txt)
# Load modules
ml Java/17.0.6 SAMtools/1.18-GCC-12.3.0
ml picard/3.3.0-Java-17
java -jar $EBROOTPICARD/picard.jar

# Global variables
PICARD=$EBROOTPICARD/picard.jar
#MARKDUPS="MarkDuplicates"
ALIGNEDFOLDER="06_bam_files"
ALIGNEDFOLDER_test="06_bam_files/test2"
METRICSFOLDER="99_metrics"

# Copy script to log folder
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98_log_files"

export JAVA_TOOL_OPTIONS="-Xms2g -Xmx50g "
export _JAVA_OPTIONS="-Xms2g -Xmx50g "

#Pass the sample number from the sbatch command
samp_num=$(($SLURM_ARRAY_TASK_ID +1))

# Fetch filename from the array
sample_name=$(cut -f1 02_info_files/SRR_Acc_List_ML.txt | sed -n "${samp_num}p")
file=${sample_name}.trimmed.fastq.gz.sorted.bam ###again need to make sure that this is a file that exists###

#removing duplicates
samtools view -f 0x2 -b $ALIGNEDFOLDER_test/$file > $ALIGNEDFOLDER_test/${sample_name}.trimmed.fastq.gz.sorted.depaired.bam
file2=${sample_name}.trimmed.fastq.gz.sorted.depaired.bam
#this will remove the duplicates stored within the same file and remove the weird reads.
#copied files over to a test folder to make sure it works first

echo "Validating sample $file2"

java -jar $PICARD ValidateSamFile \
      -I $ALIGNEDFOLDER_test/$file2 \
      -MODE SUMMARY

echo "DEduplicatING sample $file2"

java -jar $PICARD MarkDuplicates \
    -I $ALIGNEDFOLDER_test/$file2 \
    -O $ALIGNEDFOLDER_test/${sample_name}.dedup.bam \
    -METRICS_FILE $METRICSFOLDER/${sample_name}_DUP_metrics.txt \
    -VALIDATION_STRINGENCY SILENT \
    -ASSUME_SORT_ORDER queryname \
    -REMOVE_DUPLICATES true
#MarkDuplicates -I 06_bam_files/SRR19221339.trimmed.fastq.gz.sorted.bam -O 06_bam_files/SRR19221339.dedup.bam \
#-METRICS_FILE 99_metrics/SRR19221339_DUP_metrics.txt -VALIDATION_STRINGENCY SILENT -REMOVE_DUPLICATES true

echo " >>> Cleaning a bit...
"
echo "DONE! Go check your files."
