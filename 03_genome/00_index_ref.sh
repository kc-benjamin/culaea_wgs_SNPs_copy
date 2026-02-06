#!/bin/bash
#SBATCH --partition=batch
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=00_index_ref_%j.out
#SBATCH --job-name=00_index_ref 
#SBATCH --error=00_index_ref_%j.err

cd /scratch/kcb95328/Mee-Culaea-WGS

#module load nixpkgs/16.09 gcc/7.3.0

echo "Starting bwa index at: `date`" 

module load BWA/0.7.18-GCCcore-13.3.0 SAMtools/1.21-GCC-13.3.0

bwa index brook_genome_hap1_v1.fna

echo "Finished bwa index at: `date`" 

echo "Starting samtools faidx at: `date`" 

#module load samtools/1.9

samtools faidx brook_genome_hap1_v1.fna

echo "Finished samtools faidx at: `date`" 

echo "Starting picard CreateSequenceDictionary at: `date`" 

module purge
module load picard/3.3.0-Java-17
 
java -Xmx32G -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
	R=brook_genome_hap1_v1.fna \
	O=brook_genome_hap1_v1.dict

module unload nixpkgs/16.09 gcc/7.3.0 picard/2.18.9

echo "Program finished with exit code $? at: `date`" 
