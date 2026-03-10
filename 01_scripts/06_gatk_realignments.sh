#!/bin/bash

#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00-24:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out
#SBATCH --error=98_log_files/%x_%j.err
#SBATCH --array=1-1

# Copy script to log folder
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98_log_files"
cp "$SCRIPT" "$LOG_FOLDER"/"$TIMESTAMP"_"$NAME"

# CONDA_BASE=$(conda info --base)
# source ${CONDA_BASE}/etc/profile.d/conda.sh 
# conda activate /home/kcb95328/conda/envs/culaea_pkgs

# Load needed modules
module load SAMtools/1.18-GCC-12.3.0 GATK/3.8-1-Java-1.8.0_241 picard/2.18.4-Java-1.8.0_241

###conda doesnt like to work for this script since GATK is Java

# Global variables
BAM="06_bam_files"
GENOMEFOLDER="03_genome"
GENOME=$(ls -1 $GENOMEFOLDER/GCF_949316345.1_Punpun_genome.fa | xargs -n 1 basename) #changed to brook genome
GENOME_FULL="$GENOMEFOLDER/$GENOME"

# Build Bam Index
echo " >>> Realigning...
"

#Pass the sample number from the sbatch command
samp_num=$SLURM_ARRAY_TASK_ID

# Fetch filename from the array
# sample_name=$(cut -f1 02_info_files/SRR_Acc_List_ML.txt | sed -n "${samp_num}p")
sample_name=SRR19221290

#file=${sample_name}.dedup.bam ##changed to the trimmed sorted file because i dedup-ed using fastp
file=SRR19221290.dedup.bam 

echo "
     >>> Realigning TARGET for $file <<<
     "

# Index the bam file First #what is this for? to make coordinates to point the script where in the genome each read belongs
samtools sort $BAM/$file -o $BAM/${sample_name}.dedup.sorted.bam
file2=${sample_name}.dedup.sorted.bam
samtools index $BAM/$file2
#also need the bam file sorted

# Now load modules
#module purge
#module load nixpkgs/16.09 #ensure that this is right#
### THIS ONLY NEEDS TO HAPPEN ONCE####
### DO IN QLOGIN ###
# java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
#     R=$GENOMEFOLDER/$GENOME \
#     O=$GENOMEFOLDER/$GENOME.dict


# Realign
java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R $GENOMEFOLDER/$GENOME \
    -I $BAM/$file2 \
    -o $BAM/${sample_name}.intervals

echo "
     >>> Realigning INDELs for $file2 <<<
     "
java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R $GENOMEFOLDER/$GENOME \
    -I $BAM/$file2 \
    -targetIntervals $BAM/${sample_name}.intervals \
    --consensusDeterminationModel USE_READS  \
    -o $BAM/${sample_name}.realigned.bam

echo ">>> Cleaning a bit...
"
echo "
DONE! Check your files"
#conda deactivate