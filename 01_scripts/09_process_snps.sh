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

echo "Uncompressing vcf.gz file at: `date`"

gunzip -c muir_snps_filtered.vcf.gz > muir_snps_filtered_LG.vcf

awk '{gsub(/OX438541.1/,"1"); print}' muir_snps_filtered_LG.vcf > muir_snps_filtered_temp1.vcf
awk '{gsub(/OX438542.1/,"2"); print}' muir_snps_filtered_temp1.vcf > muir_snps_filtered_temp2.vcf
awk '{gsub(/OX438543.1/,"3"); print}' muir_snps_filtered_temp2.vcf > muir_snps_filtered_temp3.vcf
awk '{gsub(/OX438544.1/,"4"); print}' muir_snps_filtered_temp3.vcf > muir_snps_filtered_temp4.vcf
awk '{gsub(/OX438545.1/,"5"); print}' muir_snps_filtered_temp4.vcf > muir_snps_filtered_temp5.vcf
awk '{gsub(/OX438546.1/,"6"); print}' muir_snps_filtered_temp5.vcf > muir_snps_filtered_temp6.vcf
awk '{gsub(/OX438547.1/,"7"); print}' muir_snps_filtered_temp6.vcf > muir_snps_filtered_temp7.vcf
awk '{gsub(/OX438548.1/,"8"); print}' muir_snps_filtered_temp7.vcf > muir_snps_filtered_temp8.vcf
awk '{gsub(/OX438549.1/,"9"); print}' muir_snps_filtered_temp8.vcf > muir_snps_filtered_temp9.vcf
awk '{gsub(/OX438550.1/,"10"); print}' muir_snps_filtered_temp9.vcf > muir_snps_filtered_temp10.vcf
awk '{gsub(/OX438551.1/,"11"); print}' muir_snps_filtered_temp10.vcf > muir_snps_filtered_temp11.vcf
awk '{gsub(/OX438552.1/,"12"); print}' muir_snps_filtered_temp11.vcf > muir_snps_filtered_temp12.vcf
awk '{gsub(/OX438553.1/,"13"); print}' muir_snps_filtered_temp12.vcf > muir_snps_filtered_temp13.vcf
awk '{gsub(/OX438554.1/,"14"); print}' muir_snps_filtered_temp13.vcf > muir_snps_filtered_temp14.vcf
awk '{gsub(/OX438555.1/,"15"); print}' muir_snps_filtered_temp14.vcf > muir_snps_filtered_temp15.vcf
awk '{gsub(/OX438556.1/,"16"); print}' muir_snps_filtered_temp15.vcf > muir_snps_filtered_temp16.vcf
awk '{gsub(/OX438557.1/,"17"); print}' muir_snps_filtered_temp16.vcf > muir_snps_filtered_temp17.vcf
awk '{gsub(/OX438558.1/,"18"); print}' muir_snps_filtered_temp17.vcf > muir_snps_filtered_temp18.vcf
awk '{gsub(/OX438559.1/,"19"); print}' muir_snps_filtered_temp18.vcf > muir_snps_filtered_temp19.vcf
awk '{gsub(/OX438560.1/,"20"); print}' muir_snps_filtered_temp19.vcf > muir_snps_filtered_temp20.vcf
awk '{gsub(/OX438561.1/,"21"); print}' muir_snps_filtered_temp20.vcf > muir_snps_filtered_temp21.vcf
awk '{gsub(/OX438562.1/,"M"); print}' muir_snps_filtered_temp21.vcf > muir_snps_filtered_tempM.vcf

cp muir_snps_filtered_tempM.vcf muir_snps_filtered.vcf

rm muir_snps_filtered_temp*

echo "Converting to PLINK at: `date`"

## load all modules at once ###
#module load VCFtools/0.1.16-GCC-13.3.0 plink/1.9b_6.21-x86_64 gcta/1.26.0
#conda init bash
#conda activate /home/kcb95328/conda/envs/culaea_pkgs #this doesnt work just copy over the packages it loads
###############
# Initialize Conda (adjust the path to your conda installation if needed)
#source activate /home/kcb95328/conda/envs/culaea_pkgs
# Activate the environment
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh 
conda activate /home/kcb95328/conda/envs/culaea_pkgs

#make a chromosome map
bcftools view -H muir_snps_filtered.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > muir_snps_filtered.chrom-map.txt

#make a ped file using this chrom map
vcftools --vcf muir_snps_filtered.vcf --plink --chrom-map muir_snps_filtered.chrom-map.txt --out MU_snps

#vcftools --vcf muir_snps_filtered.vcf --plink --out MU_snps

echo "Starting PLINK filtering and removing missing data at: `date`"

#module plink/1.9b_6.21-x86_64

echo "Creating VCF file with appropriate LG labels at: `date`"

plink --file MU_snps --allow-extra-chr --recode vcf --out MU_snps

echo "Filtering at: `date`"

plink --file MU_snps --allow-extra-chr --geno 0.2 --maf 0.01 --recode --out MU_snps_geno20_maf01
#geno: removes SNPs with more than 20% missing data
#maf: removes SNPs with minor allele frequency less than 0.01
#recode: outputs in ped/map format, which is needed for GCTA

echo "Finished removing missing data at: `date`"

plink --file MU_snps_geno20_maf01 --allow-extra-chr --indep 50 5 2
plink --file MU_snps_geno20_maf01 --allow-extra-chr --extract plink.prune.in --make-bed --out MU_snps_geno20_maf01_pruned

echo "Finished pruning SNPs for LD at: `date`"

plink --file MU_snps_geno20_maf01 --allow-extra-chr --make-bed --out MU_snps_geno20_maf01

echo "Finished saving binary files for unpruned SNPs at: `date`"

#module unload plink/1.9b_6.21-x86_64

echo "Making GRM from pruned data at: `date`"

#module load gcta/1.26.0

gcta64 --bfile MU_snps_geno20_maf01_pruned --autosome --make-grm --out MU_snps_geno20_maf01_pruned

#module unload gcta/1.26.0

echo "Program finished with exit code $? at: `date`"

conda deactivate
