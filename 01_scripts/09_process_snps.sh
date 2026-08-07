#!/bin/bash
#SBATCH --partition=batch
#SBATCH --partition=batch
#SBATCH --job-name="09_process_snps"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00
#SBATCH --mem=18G
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=09_process_snps_%j.out
#SBATCH --error=09_process_snps_%j.err

echo "Converting to PLINK at: `date`"

## load all modules at once ###
#module load VCFtools/0.1.16-GCC-13.3.0 plink/1.9b_6.21-x86_64 gcta/1.26.0
###############
# Alternatively, initialize conda (adjust the path to your conda installation if needed)
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh 
conda activate /home/kcb95328/conda/envs/culaea_pkgs

#make a chromosome map
#bcftools view -H muir_snps_filtered.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > muir_snps_filtered.chrom-map.txt

#make a ped file using this chrom map
vcftools --gzvcf Ast_snps_IncFilt_filtered_renamed.vcf.gz --plink --out Ast_snps_IncFilt_filtered_renamed

#vcftools --vcf muir_snps_filtered.vcf --plink --out MU_snps
plink --file Ast_snps_IncFilt_filtered_renamed --allow-extra-chr --recode vcf --out Ast_snps_IncFilt_filtered_renamed

echo "Starting PLINK filtering and removing missing data at: `date`"

#module plink/1.9b_6.21-x86_64

echo "Creating VCF file with appropriate LG labels at: `date`"

#plink --file MU_snps --allow-extra-chr --recode vcf --out MU_snps

echo "Filtering at: `date`"

plink --file Ast_snps_IncFilt_filtered_renamed --allow-extra-chr --geno 0.2 --maf 0.01 --recode --out Ast_snps_IncFilt_filtered_renamed_geno20_maf01

#geno: removes SNPs with more than 20% missing data
#maf: removes SNPs with minor allele frequency less than 0.01
#recode: outputs in ped/map format, which is needed for GCTA

#5/18/2026: I wanted to see what kind of plot is made by plink at this stage; see here:
plink --file Muir_snps_geno20_maf01 --assoc perm fisher --allow-extra-chr --allow-no-sex --pheno SL_pheno_numbers.txt --out Muir_snps_geno20_maf01_assoc

### 5/26/2026: trying this on the amhy- males only
plink --file Shunda_snps_geno20_maf01 --assoc perm fisher --allow-extra-chr --allow-no-sex --pheno SL_geno_SRmales_and_females_only.txt --out Shunda_amhy-males_and_females

#### 5/27/2026: for Amhy Masked Shunda, all males and females
plink --file Shunda_AmhyMasked_snps_geno20_maf01 --assoc perm fisher --allow-extra-chr --allow-no-sex --pheno SL_pheno_numbers.txt --out Shunda_AmhyMasked_snps_geno20_maf01_assoc

### 5/29/2026: for Muir
plink --file Muir_snps_geno20_maf01 --assoc perm fisher --allow-extra-chr --allow-no-sex --pheno ML_sex_pheno_numbers.txt --out Muir_snps_geno20_maf01_assoc

### 6/1/2026: Muir
# for amhy- males and amhy- females
plink --file Muir_snps_geno20_maf01 --assoc perm fisher --allow-extra-chr --allow-no-sex --pheno ML_all_phenotypes.txt --pheno-name Amhy_Null --out Muir_snps_amhy_null

plink --file Astotin_snps_geno20_maf01 --assoc perm fisher --allow-extra-chr --allow-no-sex --pheno AS_sex_pheno_amhynull-only.txt --out Astotin_snps_amhy_null

### 6/3/2026: Muir but with all chromosomes
plink --file Muir_snps_ALL_geno20_maf01 --assoc perm fisher --allow-extra-chr --allow-no-sex --pheno ML_sex_pheno_numbers.txt --out Muir_snps_ALL
plink --file Muir_snps_ALL_geno20_maf01 --assoc perm fisher --allow-extra-chr --allow-no-sex --pheno ML_all_phenotypes.txt --pheno-name Amhy_Null --out Muir_snps_ALL_amhy_null

### 6/9/2026: on Astotin but with increased filtering at the samtools step
plink --file Ast_snps_IncFilt_filtered_renamed_geno20_maf01 --assoc perm fisher --allow-extra-chr --allow-no-sex --pheno AS_sex_pheno_numbers.txt --out Ast_snps_IncFilt_MFall
plink --file Ast_snps_IncFilt_filtered_renamed_geno20_maf01 --assoc perm fisher --allow-extra-chr --allow-no-sex --pheno AS_sex_pheno_amhynull-only.txt --out Ast_snps_IncFilt_amhynull




echo "Finished removing missing data at: `date`"

plink --file MU_snps_geno20_maf01 --allow-extra-chr --indep 50 5 2
plink --file MU_snps_geno20_maf01 --make-pheno test-pheno.txt 2 --mpheno 4 --update-sex test-sex.txt --allow-extra-chr --extract plink.prune.in --make-bed --out MU_snps_geno20_maf01_pruned

echo "Finished pruning SNPs for LD at: `date`"

plink --file MU_snps_geno20_maf01 --update-chr brook_chr_concat.txt --make-pheno test-pheno.txt 2 --mpheno 4 --update-sex test-sex.txt --allow-extra-chr --make-bed --out MU_snps_geno20_maf01
#plink --bfile MU_snps_geno20_maf01 --make-pheno test-pheno.txt 2 --mpheno 4 --update-sex test-sex.txt --allow-extra-chr --out MU_snps_geno20_maf01_pheno

echo "Finished saving binary files for unpruned SNPs at: `date`"

#make a phenotype file: do this using qlogin, then continue
#echo "Files are ready. Use qlogin to make a phenotype file and then make the GRM using the following steps in script 09"

#plink --pheno 02_info_files/SraRunTable-metadata-ML_sexincl.csv --allow-extra-chr --pheno-name "sex_genotype"

#module unload plink/1.9b_6.21-x86_64

echo "Making GRM from pruned data at: `date`"

#module load gcta/1.26.0

#gcta64 --bfile MU_snps_geno20_maf01_pruned --autosome --make-grm --out MU_snps_geno20_maf01_pruned

#module unload gcta/1.26.0

#echo "Program finished with exit code $? at: `date`"

conda deactivate
