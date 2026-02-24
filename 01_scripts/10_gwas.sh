#!/bin/bash

#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=1024M
#SBATCH --time=0-12:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out
#SBATCH --error=98_log_files/%x_%j.err

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh 
conda activate /home/kcb95328/conda/envs/culaea_pkgs

#genome-wide complex trait analysis (GCTA) using the mlma (mixed linear model association) method
gcta64 --mlma --bfile MU_snps_geno20_maf01 --grm MU_snps_geno20_maf01_pruned --pheno MU_snps.phen --out MU_brook_geno20_maf01

conda deactivate