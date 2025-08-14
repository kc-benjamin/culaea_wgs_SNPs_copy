#!/bin/bash
#SBATCH --partition=batch
#SBATCH --time=06:00:00
#SBATCH --mem=30G
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=09_process_snps_%j.out
#SBATCH --job-name=09_process_snps 

echo "Uncompressing vcf.gz file at: `date`"

gunzip -c shunda_snps_filtered.vcf.gz > shunda_snps_filtered_LG.vcf

awk '{gsub(/OX438541.1/,"1"); print}' shunda_snps_filtered_LG.vcf > shunda_snps_filtered_temp1.vcf
awk '{gsub(/OX438542.1/,"2"); print}' shunda_snps_filtered_temp1.vcf > shunda_snps_filtered_temp2.vcf
awk '{gsub(/OX438543.1/,"3"); print}' shunda_snps_filtered_temp2.vcf > shunda_snps_filtered_temp3.vcf
awk '{gsub(/OX438544.1/,"4"); print}' shunda_snps_filtered_temp3.vcf > shunda_snps_filtered_temp4.vcf
awk '{gsub(/OX438545.1/,"5"); print}' shunda_snps_filtered_temp4.vcf > shunda_snps_filtered_temp5.vcf
awk '{gsub(/OX438546.1/,"6"); print}' shunda_snps_filtered_temp5.vcf > shunda_snps_filtered_temp6.vcf
awk '{gsub(/OX438547.1/,"7"); print}' shunda_snps_filtered_temp6.vcf > shunda_snps_filtered_temp7.vcf
awk '{gsub(/OX438548.1/,"8"); print}' shunda_snps_filtered_temp7.vcf > shunda_snps_filtered_temp8.vcf
awk '{gsub(/OX438549.1/,"9"); print}' shunda_snps_filtered_temp8.vcf > shunda_snps_filtered_temp9.vcf
awk '{gsub(/OX438550.1/,"10"); print}' shunda_snps_filtered_temp9.vcf > shunda_snps_filtered_temp10.vcf
awk '{gsub(/OX438551.1/,"11"); print}' shunda_snps_filtered_temp10.vcf > shunda_snps_filtered_temp11.vcf
awk '{gsub(/OX438552.1/,"12"); print}' shunda_snps_filtered_temp11.vcf > shunda_snps_filtered_temp12.vcf
awk '{gsub(/OX438553.1/,"13"); print}' shunda_snps_filtered_temp12.vcf > shunda_snps_filtered_temp13.vcf
awk '{gsub(/OX438554.1/,"14"); print}' shunda_snps_filtered_temp13.vcf > shunda_snps_filtered_temp14.vcf
awk '{gsub(/OX438555.1/,"15"); print}' shunda_snps_filtered_temp14.vcf > shunda_snps_filtered_temp15.vcf
awk '{gsub(/OX438556.1/,"16"); print}' shunda_snps_filtered_temp15.vcf > shunda_snps_filtered_temp16.vcf
awk '{gsub(/OX438557.1/,"17"); print}' shunda_snps_filtered_temp16.vcf > shunda_snps_filtered_temp17.vcf
awk '{gsub(/OX438558.1/,"18"); print}' shunda_snps_filtered_temp17.vcf > shunda_snps_filtered_temp18.vcf
awk '{gsub(/OX438559.1/,"19"); print}' shunda_snps_filtered_temp18.vcf > shunda_snps_filtered_temp19.vcf
awk '{gsub(/OX438560.1/,"20"); print}' shunda_snps_filtered_temp19.vcf > shunda_snps_filtered_temp20.vcf
awk '{gsub(/OX438561.1/,"21"); print}' shunda_snps_filtered_temp20.vcf > shunda_snps_filtered_temp21.vcf
awk '{gsub(/OX438562.1/,"M"); print}' shunda_snps_filtered_temp21.vcf > shunda_snps_filtered_tempM.vcf

cp shunda_snps_filtered_tempM.vcf shunda_snps_filtered.vcf

rm shunda_snps_filtered_temp*

echo "Converting to PLINK at: `date`"

module vcftools/0.1.16

vcftools --vcf shunda_snps_filtered.vcf --plink --out SHU_snps

echo "Starting PLINK filtering and removing missing data at: `date`"

module plink/1.9b_6.21-x86_64

echo "Creating VCF file with appropriate LG labels at: `date`"

plink --file SHU_snps --allow-extra-chr --recode vcf --out SHU_snps

echo "Filtering at: `date`"

plink --file SHU_snps --geno 0.2 --maf 0.01 --recode --out SHU_snps_geno20_maf01

echo "Finished removing missing data at: `date`"

plink --file SHU_snps_geno20_maf01 --indep 50 5 2
plink --file SHU_snps_geno20_maf01 --extract plink.prune.in --make-bed --out SHU_snps_geno20_maf01_pruned

echo "Finished pruning SNPs for LD at: `date`"

plink --file SHU_snps_geno20_maf01 --make-bed --out SHU_snps_geno20_maf01

echo "Finished saving binary files for unpruned SNPs at: `date`"

module unload plink/1.9b_6.21-x86_64

echo "Making GRM from pruned data at: `date`"

module load gcta/1.26.0

gcta64 --bfile SHU_snps_geno20_maf01_pruned --autosome --make-grm --out SHU_snps_geno20_maf01_pruned

module unload gcta/1.26.0

echo "Program finished with exit code $? at: `date`"
