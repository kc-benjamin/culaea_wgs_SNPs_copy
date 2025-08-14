#!/bin/bash
#SBATCH --partition=batch
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=ref_pull_ncbi_%j.out
#SBATCH --job-name=00_index_ref 
#SBATCH --error=ref_pull_ncbi_%j.err

module load datasets

datasets download genome accession PRJEB33823 --filename punpun-ref-genome.fa