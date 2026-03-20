#!/bin/bash
#SBATCH --partition=batch
#SBATCH --partition=batch
#SBATCH --job-name="10a_shunda_plink"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=2G
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=10a_shunda_plink_%j.out

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh 
conda activate /home/kcb95328/conda/envs/culaea_pkgs

echo "step 1: convert vcf to plink format"
vcftools --vcf shunda_snps_renamed.vcf --plink --out SHU_snps2

echo "step 2: convert plink to vcf format"
plink --file SHU_snps2 --recode vcf --out SHU_snps2

echo "step 3: filter plink file for missingness and minor allele frequency"
plink --file SHU_snps2 --geno 0.2 --maf 0.01 --recode --out SHU_snps_geno20_maf01_2

echo "step 4: prune for linkage disequilibrium"
plink --file SHU_snps_geno20_maf01_2 --indep 50 5 2 

echo "step 5: create pruned plink file"
plink --file SHU_snps_geno20_maf01_2 --extract plink.prune.in --make-bed --out SHU_snps_geno20_maf01_2_pruned

echo "step 6: make bed file for filtered file"
plink --file SHU_snps_geno20_maf01_2 --make-bed --out SHU_snps_geno20_maf01_2

echo "step 7: create genetic relationship matrix"
gcta64 --bfile SHU_snps_geno20_maf01_2_pruned --make-grm --out SHU_snps_geno20_maf01_2_pruned2

echo "step 8: run the mlma"
gcta64 --mlma --bfile SHU_snps_geno20_maf01_2_pruned --grm SHU_snps_geno20_maf01_2_pruned --pheno SHU_snps.phen --out SHU_snps_jon_phenos_spines_2
