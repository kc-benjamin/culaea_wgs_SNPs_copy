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


module load GCTA/1.94.1-gfbf-2023a

gcta64 --bfile MU_snps_geno20_maf01_pruned --autosome-num 22 --autosome --make-grm --out MU_snps_geno20_maf01_pruned

gcta64 --mlma --bfile MU_snps_geno20_maf01 --grm MU_snps_geno20_maf01_pruned --pheno MU_snps.phen --out MU_brook_geno20_maf01

module unload GCTA/1.94.1-gfbf-2023a