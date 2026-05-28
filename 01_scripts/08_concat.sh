#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name="08_concat"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00-12:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out
#SBATCH --error=98_log_files/%x_%j.err

module load VCFtools/0.1.16-GCC-13.3.0 BCFtools/1.21-GCC-13.3.0 tabix/0.2.6-GCCcore-13.3.0 
#module intel/2020.1.217 
OUT="/scratch/kcb95328/MuirLakeBrooks/08_PLINK_new"

vcf-concat $(ls -1 /scratch/kcb95328/MuirLakeBrooks/07_vcfs_new/* | perl -pe 's/\n/ /g') > $OUT/Muir_snps.vcf

bgzip $OUT/Muir_snps.vcf
tabix --force --preset vcf $OUT/Muir_snps.vcf

bcftools filter -e 'MQ < 30' $OUT/Muir_snps.vcf.gz -Oz > $OUT/tmp.vcf.gz

vcftools --gzvcf $OUT/tmp.vcf.gz --minQ 30 --minGQ 20 --minDP 1 --recode --recode-INFO-all --stdout > $OUT/Muir_snps_filtered.vcf

rm tmp.vcf.gz

bgzip $OUT/Muir_snps_filtered.vcf
tabix -p vcf $OUT/Muir_snps_filtered.vcf.gz