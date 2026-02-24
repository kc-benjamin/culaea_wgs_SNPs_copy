#!/bin/bash
#SBATCH --partition=batch
#SBATCH --partition=batch
#SBATCH --job-name="09a_gcta"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00
#SBATCH --mem=10G
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=09a_gcta_%j.out


module load gcta/1.26.0

gcta64 --bfile MU_snps_geno20_maf01_pruned --autosome-num 22 --autosome --make-grm --out MU_snps_geno20_maf01_pruned

module unload gcta/1.26.0