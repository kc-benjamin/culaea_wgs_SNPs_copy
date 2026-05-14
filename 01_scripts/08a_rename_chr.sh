#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name="08a_rename_chr"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00-24:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out
#SBATCH --error=98_log_files/%x_%j.err

ml BCFtools/1.21-GCC-13.3.0 tabix/0.2.6-GCCcore-13.3.0
gunzip -c Astotin_snps_filtered.vcf.gz > Astotin_snps_filtered.vcf

awk '{gsub(/PGA_scaffold22__123_contigs__length_35631191/,"4"); print}' Astotin_snps_filtered.vcf > Astotin_snps_filtered_temp1.vcf
awk '{gsub(/PGA_scaffold23__90_contigs__length_32925445/,"13"); print}' Astotin_snps_filtered_temp1.vcf > Astotin_snps_filtered_temp2.vcf
awk '{gsub(/PGA_scaffold21__64_contigs__length_27087477/,"2"); print}' Astotin_snps_filtered_temp2.vcf > Astotin_snps_filtered_temp3.vcf
awk '{gsub(/PGA_scaffold20__81_contigs__length_24420605/,"9"); print}' Astotin_snps_filtered_temp3.vcf > Astotin_snps_filtered_temp4.vcf
awk '{gsub(/PGA_scaffold18__88_contigs__length_23069817/,"12"); print}' Astotin_snps_filtered_temp4.vcf > Astotin_snps_filtered_temp5.vcf
awk '{gsub(/PGA_scaffold19__71_contigs__length_22530169/,"17"); print}' Astotin_snps_filtered_temp5.vcf > Astotin_snps_filtered_temp6.vcf
awk '{gsub(/PGA_scaffold17__83_contigs__length_22236360/,"8"); print}' Astotin_snps_filtered_temp6.vcf > Astotin_snps_filtered_temp7.vcf
awk '{gsub(/PGA_scaffold14__88_contigs__length_21401847/,"20"); print}' Astotin_snps_filtered_temp7.vcf > Astotin_snps_filtered_temp8.vcf
awk '{gsub(/PGA_scaffold16__56_contigs__length_21394868/,"1"); print}' Astotin_snps_filtered_temp8.vcf > Astotin_snps_filtered_temp9.vcf
awk '{gsub(/PGA_scaffold15__61_contigs__length_21329229/,"19"); print}' Astotin_snps_filtered_temp9.vcf > Astotin_snps_filtered_temp10.vcf
awk '{gsub(/PGA_scaffold12__56_contigs__length_20552211/,"6"); print}' Astotin_snps_filtered_temp10.vcf > Astotin_snps_filtered_temp11.vcf
awk '{gsub(/PGA_scaffold10__89_contigs__length_20417627/,"16"); print}' Astotin_snps_filtered_temp11.vcf > Astotin_snps_filtered_temp12.vcf
awk '{gsub(/PGA_scaffold13__42_contigs__length_19371112/,"7"); print}' Astotin_snps_filtered_temp12.vcf > Astotin_snps_filtered_temp13.vcf
awk '{gsub(/PGA_scaffold11__71_contigs__length_19041060/,"11"); print}' Astotin_snps_filtered_temp13.vcf > Astotin_snps_filtered_temp14.vcf
awk '{gsub(/PGA_scaffold7__70_contigs__length_18533579/,"15"); print}' Astotin_snps_filtered_temp14.vcf > Astotin_snps_filtered_temp15.vcf
awk '{gsub(/PGA_scaffold9__59_contigs__length_18519828/,"21"); print}' Astotin_snps_filtered_temp15.vcf > Astotin_snps_filtered_temp16.vcf
awk '{gsub(/PGA_scaffold8__67_contigs__length_18086311/,"10"); print}' Astotin_snps_filtered_temp16.vcf > Astotin_snps_filtered_temp17.vcf
awk '{gsub(/PGA_scaffold5__81_contigs__length_16821318/,"14"); print}' Astotin_snps_filtered_temp17.vcf > Astotin_snps_filtered_temp18.vcf
awk '{gsub(/PGA_scaffold6__72_contigs__length_16467905/,"5"); print}' Astotin_snps_filtered_temp18.vcf > Astotin_snps_filtered_temp19.vcf
awk '{gsub(/PGA_scaffold4__53_contigs__length_16101261/,"3"); print}' Astotin_snps_filtered_temp19.vcf > Astotin_snps_filtered_temp20.vcf
awk '{gsub(/PGA_scaffold3__61_contigs__length_14436489/,"7"); print}' Astotin_snps_filtered_temp20.vcf > Astotin_snps_filtered_temp21.vcf
awk '{gsub(/PGA_scaffold2__66_contigs__length_13522551/,"18"); print}' Astotin_snps_filtered_temp21.vcf > Astotin_snps_filtered_temp22.vcf
awk '{gsub(/PGA_scaffold1__43_contigs__length_10481396/,"1"); print}' Astotin_snps_filtered_temp22.vcf > Astotin_snps_filtered_renamed.vcf

rm Astotin_snps_filtered_temp*.vcf

bgzip Astotin_snps_filtered_renamed.vcf
tabix -p vcf Astotin_snps_filtered_renamed.vcf.gz

