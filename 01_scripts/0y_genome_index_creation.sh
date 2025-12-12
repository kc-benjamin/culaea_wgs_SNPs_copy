#!/bin/bash

#SBATCH --partition=batch  
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=0-24:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/genome_index_%j.out

module load BWA/0.7.18-GCCcore-13.3.0 SAMtools/1.18-GCC-12.3.0

#location variables
GENOMEFOLDER="03_genome"
GENOME=$(ls -1 $GENOMEFOLDER/brook_genome_hap1_v1.fa | xargs -n 1 basename) #changed to brook genome
GENOME_FULL="$GENOMEFOLDER/$GENOME"

#bwa index -p brook_genome_hap1_v1.fa -a is $GENOME_FULL

samtools faidx $GENOME_FULL