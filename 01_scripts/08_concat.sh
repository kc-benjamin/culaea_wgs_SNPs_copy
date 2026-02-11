#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name="08_concat"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00-23:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out
#SBATCH --error=98_log_files/%x_%j.err

module load VCFtools/0.1.16-GCC-13.3.0 BCFtools/1.21-GCC-13.3.0 tabix/0.2.6-GCCcore-13.3.0
#module intel/2020.1.217 

vcf-concat $(ls -1 ./07_raw_VCFs/output.vcf | perl -pe 's/\n/ /g') > muir_snps.vcf

bgzip muir_snps.vcf
tabix -p vcf muir_snps.vcf.gz

bcftools filter -e 'MQ < 30' muir_snps.vcf.gz -Oz > tmp.vcf.gz

vcftools --gzvcf tmp.vcf.gz --minQ 30 --minGQ 20 --minDP 5 --max-alleles 2 --recode --recode-INFO-all --stdout > muir_snps_filtered.vcf

rm tmp.vcf.gz

bgzip muir_snps_filtered.vcf
tabix -p vcf muir_snps_filtered.vcf.gz
