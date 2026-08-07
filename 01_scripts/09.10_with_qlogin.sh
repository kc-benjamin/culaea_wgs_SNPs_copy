#cleaned up plink and gcta script
#designed for use with qlogin

#1) activate conda
conda activate /home/kcb95328/conda/envs/culaea_pkgs

#2) make a chromosome map
bcftools view -H muir_snps_filtered.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > muir_snps_filtered.chrom-map.txt

#3) make a ped file using this chrom map
vcftools --vcf muir_snps_filtered.vcf --plink --chrom-map muir_snps_filtered.chrom-map.txt --out MU_snps

#4) recode for plink
plink --file MUI_snps --allow-extra-chr --recode-vcf --out MUI_snps

#5) geno and maf filtering
plink --file MU_snps --pheno MUI_snps.phen --allow-extra-chr --geno 0.2 --maf 0.01 --recode --out MU_snps_geno20_maf01

#6) prune snps
plink --file MUI_snps_geno20_maf01 --indep 50 5 2
plink --file MUI_snps_geno20_maf01 --extract plink.prune.in --make-bed --out MUI_snps_geno20_maf01_pruned

#7) make bed files for the step prior
plink --file MUI_snps_geno20_maf01 --make-bed --out MUI_snps_geno20_maf01

#8) make the grm
gcta64 --bfile MUI_snps_geno20_maf01_pruned --autosome --make-grm --out MUI_snps_geno20_maf01_pruned

#9) run the mlma
gcta64 --mlma --bfile SU_snps_geno20_maf01_pruned --grm SU_snps_geno20_maf01_pruned --pheno SHU_snps.phen --out SHU_snps_jon_phenos_spines



