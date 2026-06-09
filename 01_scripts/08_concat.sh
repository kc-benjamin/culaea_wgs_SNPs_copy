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
DIR="/scratch/kcb95328/IncreasedFiltering-Astotin"
VCFS="/07_vcfs"
PLINK="/08_PLINK"

vcf-concat $(ls -1 ${DIR}/${VCFS}/* | perl -pe 's/\n/ /g') > $DIR/$PLINK/Ast_snps_IncFilt.vcf

bgzip $DIR/$PLINK/Ast_snps_IncFilt.vcf
tabix --force --preset vcf $DIR/$PLINK/Ast_snps_IncFilt.vcf

bcftools filter -e 'MQ < 30' $DIR/$PLINK/Ast_snps_IncFilt.vcf.gz -Oz > $DIR/$PLINK/tmp.vcf.gz
#mapping quality is 30 here--bams only align in the "IncreasedFiltering" folders if MQ > 20

vcftools --gzvcf $DIR/$PLINK/tmp.vcf.gz --minQ 30 --minGQ 30 --minDP 1 --recode --recode-INFO-all --stdout > $DIR/$PLINK/Ast_snps_IncFilt_filtered.vcf
#from Jon: minimum quality is 30, minimum genotype quality is 20, and minimum depth is 1
#I increased the GQ (confidence that a genotype call is correct) to 30

rm $DIR/$PLINK/tmp.vcf.gz

bgzip $DIR/$PLINK/Ast_snps_IncFilt_filtered.vcf
tabix -p vcf $DIR/$PLINK/Ast_snps_IncFilt_filtered.vcf.gz