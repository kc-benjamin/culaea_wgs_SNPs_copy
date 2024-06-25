#!/bin/bash
#SBATCH --job-name="08_concat"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00-23:00:00
#SBATCH --account=def-jonmee
#SBATCH --mail-user=jmee@mtroyal.ca
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out

module load vcftools bcftools
module load StdEnv/2020 intel/2020.1.217 tabix/0.2.6

vcf-concat $(ls -1 ./07_raw_VCFs/* | perl -pe 's/\n/ /g') > shunda_snps.vcf

bgzip shunda_snps.vcf
tabix -p vcf shunda_snps.vcf.gz

bcftools filter -e 'MQ < 30' shunda_snps.vcf.gz -Oz > tmp.vcf.gz

vcftools --gzvcf tmp.vcf.gz --minQ 30 --minGQ 20 --minDP 5 --max-alleles 2 --recode --recode-INFO-all --stdout > shunda_snps_filtered.vcf

rm tmp.vcf.gz

bgzip shunda_snps_filtered.vcf
tabix -p vcf shunda_snps_filtered.vcf.gz
