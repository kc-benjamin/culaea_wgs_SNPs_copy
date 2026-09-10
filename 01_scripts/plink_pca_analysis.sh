#!/bin/bash
#SBATCH --partition=batch
#SBATCH --partition=batch
#SBATCH --job-name="plink_pca_analysis"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00
#SBATCH --mem=18G
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=plink_pca_analysis_%j.out
#SBATCH --error=plink_pca_analysis_%j.err

#writing space for PLINK analysis using covariate structure baked into commands
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate /home/kcb95328/conda/envs/culaea_pkgs

#this code assumes you have already generated the filtered, renamed, maf and genotype files

#plink --file Shunda_snps_geno20_maf01 --pca 10 --out Shunda_pca_out
plink --file Shunda_snps_geno20_maf01 --logistic mperm=10000 --allow-extra-chr --allow-no-sex --covar Shunda_pca_out.eigenvec --covar-number 1-4 --pheno /home/kcb95328/Info-Shunda/SL_pheno_numbers_all.txt --mpheno 1 --out Shunda_with_pca_GWAS_results
