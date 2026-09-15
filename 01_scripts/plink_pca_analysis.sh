#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name="plink_pca_analysis_Astotin"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=012:00:00
#SBATCH --mem=256G
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=plink_pca_analysis_%j.out
#SBATCH --error=plink_pca_analysis_%j.err

#writing space for PLINK analysis using covariate structure baked into commands
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate /home/kcb95328/conda/envs/culaea_pkgs

#this code assumes you have already generated the filtered, renamed, maf and genotype files
cd $SLURM_SUBMIT_DIR #should be the same as your input folder
INPUT="/scratch/kcb95328/AstotinLakeBrooks/"
DATA="/home/kcb95328/Info-Astotin/"


plink --file Astotin_snps_geno20_maf01 --pca 10 --out Astotin_pca_out
plink --file Astotin_snps_geno20_maf01 --logistic mperm=10000 --allow-extra-chr --allow-no-sex --covar Astotin_pca_out.eigenvec --covar-number 1-4 --pheno $DATA/AS_sex_pheno_numbers.txt --mpheno 1 --out Astotin_with_pca_GWAS_results
#plink --file Astotin_snps_geno20_maf01 --logistic mperm=10000 --allow-extra-chr --allow-no-sex --covar Astotin_pca_out.eigenvec --covar-number 1-4 --pheno $DATA/AL_pheno_numbers_all.txt --mpheno 1 --out Astotin_with_pca_GWAS_results

conda deactivate