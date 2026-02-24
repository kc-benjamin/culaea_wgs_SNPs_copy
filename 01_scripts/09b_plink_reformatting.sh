#!/bin/bash
#SBATCH --partition=batch
#SBATCH --partition=batch
#SBATCH --job-name="09b_plink_reformatting"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00
#SBATCH --mem=2G
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=09b_plink_reformatting_%j.out

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh 
conda activate /home/kcb95328/conda/envs/culaea_pkgs

plink --file MU_snps_geno20_maf01 --update-chr brook_chr_for_plink.txt --make-pheno test-pheno.txt 2 --mpheno 4 --update-sex test-sex.txt --allow-extra-chr --make-bed --out MU_snps_geno20_maf01_scaffolds

conda deactivate