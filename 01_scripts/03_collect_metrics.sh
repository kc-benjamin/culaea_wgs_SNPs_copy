#!/bin/bash

#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=0-06:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out
#SBATCH --array=0-96
PREFIX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" 02_info_files/SRR_Acc_List_ML.txt)
# Load modules
module load Java/11.0.20 picard/2.25.1-Java-ll R/4.4.2-gfbf-2024a

# Global variables
GENOMEFOLDER="/scratch/kcb95328/Mee-Culaea-WGS/03_genome"
GENOME=$(ls -1 $GENOMEFOLDER/*fa | xargs -n 1 basename)
ALIGNEDFOLDER="/scratch/kcb95328/Mee-Culaea-WGS/06_bam_files"
METRICSFOLDER="/scratch/kcb95328/Mee-Culaea-WGS/99_metrics"
PICARD=$EBROOTPICARD/picard.jar
ALIGN="CollectAlignmentSummaryMetrics"
INSERT="CollectInsertSizeMetrics"
COVERAGE="CollectWgsMetricsWithNonZeroCoverage"
# BAM_FILE="02_info_files/bam_split/bam0000"

# Copy script to log folder
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
SCRIPTNAME=$(basename $0)
LOG_FOLDER="/scratch/kcb95328/Mee-Culaea-WGS/98_log_files"
cp $SCRIPT $LOG_FOLDER/${TIMESTAMP}_${SCRIPTNAME}

#Pass the sample number from the sbatch command
samp_num=PREFIX

    # Fetch filename from the array
    file=$(cut -f1 /home/kcb95328/culaea_wgs_SNPs_copy/02_info_files/SRR_Acc_List_ML.txt | sed -n "${samp_num}p")
    bamfile=${file}.sorted.bam

    echo \n">>> Computing alignment metrics for $file <<<"\n
    java -jar $PICARD $ALIGN \
        R=$GENOMEFOLDER/$GENOME \
        I=$ALIGNEDFOLDER/$bamfile \
        O=$METRICSFOLDER/${file}_alignment_metrics.txt

    echo \n">>> Computing insert size metrics for $file <<<"\n
    java -jar $PICARD $INSERT \
        I=$ALIGNEDFOLDER/$bamfile \
        OUTPUT=$METRICSFOLDER/${file}_insert_size_metrics.txt \
        HISTOGRAM_FILE=$METRICSFOLDER/${file}_insert_size_histogram.pdf

    echo \n">>> Computing coverage metrics for $file <<<"\n
    java -jar $PICARD $COVERAGE \
        R=$GENOMEFOLDER/$GENOME \
        I=$ALIGNEDFOLDER/$bamfile \
        OUTPUT=$METRICSFOLDER/${file}_collect_wgs_metrics.txt\
        CHART=$METRICSFOLDER/${file}_collect_wgs_metrics.pdf

    echo \n">>> DONE! <<<"\n
