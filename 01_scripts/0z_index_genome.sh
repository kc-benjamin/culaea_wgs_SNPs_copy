#!/bin/bash

#SBATCH --partition=batch  
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=0-24:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/0z_genome_%j.out
#SBATCH --error=98_log_files/0z_genome_%j.err

ml BWA/0.7.18-GCCcore-13.3.0
bwa index -p brook_genome_hap1_v1.fa -a is brook_genome_hap1_v1.fa