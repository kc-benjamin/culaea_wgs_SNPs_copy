##writing space for 09 plink

conda activate /home/kcb95328/conda/envs/culaea_pkgs

gunzip -c shunda_snps_filtered.vcf.gz > shunda_snps_filtered.vcf

vcftools --vcf shunda_snps_renamed.vcf --plink --out SHU_snps
#Unrecognized values used for CHROM: PGA_scaffold10__89_contigs__length_20417627 - Replacing with 0.
#interrupted after this error
vcftools --vcf shunda_snps_filtered.vcf --plink --out SHU_snps --allow-extra-chr

#need to make a chr map first
bcftools view -H shunda_snps_filtered.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > shunda_snps_filtered.chrom-map.txt
#back to the first step:
vcftools --vcf shunda_snps_filtered.vcf --plink --chrom-map shunda_snps_filtered.chrom-map.txt --out SU_snps

plink --file SHU_snps --recode vcf --out SHU_snps

plink --file SHU_snps --geno 0.2 --maf 0.01 --recode --out SU_snps_geno20_maf01
#plink --file SU_snps --geno 0.2 --maf 0.01 --recode --allow-extra-chr --out SU_snps_geno20_maf01
plink --file SHU_snps --geno 0.1 --maf 0.01 --recode --out SU_snps_geno10_maf01


plink --file SU_snps_geno20_maf01 --indep 50 5 2 
plink --file SU_snps_geno20_maf01 --indep 50 5 2 --allow-extra-chr
plink --file SU_snps_geno10_maf01 --indep 50 5 2 


plink --file SU_snps_geno20_maf01 --extract plink.prune.in --make-bed --out SU_snps_geno20_maf01_pruned

plink --file SU_snps_geno20_maf01 --make-bed --out SU_snps_geno20_maf01

gcta64 --bfile SU_snps_geno20_maf01_pruned --make-grm --out SU_snps_geno20_maf01_pruned2
## for use with Jon's grm/plink files from Dryad
gcta64 --mlma --bfile SHU_snps_geno20_maf01 --grm SHU_snps_geno20_maf01_pruned --pheno sample_ID_sex_SHU.phen --out SHU_snps_brook_jons_vcf
