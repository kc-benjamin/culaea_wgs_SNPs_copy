#!/bin/bash
#SBATCH --partition=batch
#SBATCH --partition=batch
#SBATCH --job-name="star-align"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=1-00:00:00
#SBATCH --mem=24G
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=99_log_files/star-align_%j.out
#SBATCH --error=99_log_files/star-align_%j.err

cd $SLURM_SUBMIT_DIR
ml STAR/2.7.11a-GCC-12.3.0

#variables
GENOMEDIR="02_genome"
FASTA="brook_genome_hap1_v1.fasta"
GTF="brook_genome_hap1_v1.gff"

STAR --runThreadN $SLURM_CPUS_PER_TASK \
     --runMode genomeGenerate \
     --genomeDir $GENOMEDIR \
     --genomeSAindexNbases 13 \
     --genomeFastaFiles $SLURM_SUBMIT_DIR/brook_genome_hap1_v1.fasta \
     --sjdbGTFfile $SLURM_SUBMIT_DIR/brook_genome_hap1_v1.gff \