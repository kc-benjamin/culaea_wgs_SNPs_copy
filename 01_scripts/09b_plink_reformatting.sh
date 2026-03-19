#!/bin/bash
#SBATCH --partition=batch
#SBATCH --partition=batch
#SBATCH --job-name="09b_plink_reformatting"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --mem=1G
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=09b_plink_reformatting_%j.out

awk '{gsub(/PGA_scaffold22__123_contigs__length_35631191/,"4"); print}' shunda_snps_filtered_LG.vcf > shunda_snps_filtered_temp1.vcf
awk '{gsub(/PGA_scaffold21__64_contigs__length_27087477/,"2"); print}' shunda_snps_filtered_temp1.vcf > shunda_snps_filtered_temp2.vcf
awk '{gsub(/PGA_scaffold20__81_contigs__length_24420605/,"9"); print}' shunda_snps_filtered_temp2.vcf > shunda_snps_filtered_temp3.vcf
awk '{gsub(/PGA_scaffold18__88_contigs__length_23069817/,"12"); print}' shunda_snps_filtered_temp3.vcf > shunda_snps_filtered_temp4.vcf
awk '{gsub(/PGA_scaffold19__71_contigs__length_22530169/,"17"); print}' shunda_snps_filtered_temp4.vcf > shunda_snps_filtered_temp5.vcf
awk '{gsub(/PGA_scaffold17__83_contigs__length_22236360/,"8"); print}' shunda_snps_filtered_temp5.vcf > shunda_snps_filtered_temp6.vcf
awk '{gsub(/PGA_scaffold14__88_contigs__length_21401847/,"20"); print}' shunda_snps_filtered_temp6.vcf > shunda_snps_filtered_temp7.vcf
awk '{gsub(/PGA_scaffold16__56_contigs__length_21394868/,"1"); print}' shunda_snps_filtered_temp7.vcf > shunda_snps_filtered_temp8.vcf
awk '{gsub(/PGA_scaffold15__61_contigs__length_21329229/,"19"); print}' shunda_snps_filtered_temp8.vcf > shunda_snps_filtered_temp9.vcf
awk '{gsub(/PGA_scaffold12__56_contigs__length_20552211/,"6"); print}' shunda_snps_filtered_temp9.vcf > shunda_snps_filtered_temp10.vcf
awk '{gsub(/PGA_scaffold10__89_contigs__length_20417627/,"16"); print}' shunda_snps_filtered_temp10.vcf > shunda_snps_filtered_temp11.vcf
awk '{gsub(/PGA_scaffold13__42_contigs__length_19371112/,"7"); print}' shunda_snps_filtered_temp11.vcf > shunda_snps_filtered_temp12.vcf
awk '{gsub(/PGA_scaffold11__71_contigs__length_19041060/,"11"); print}' shunda_snps_filtered_temp12.vcf > shunda_snps_filtered_temp13.vcf
awk '{gsub(/PGA_scaffold7__70_contigs__length_18533579/,"15"); print}' shunda_snps_filtered_temp13.vcf > shunda_snps_filtered_temp14.vcf
awk '{gsub(/PGA_scaffold9__59_contigs__length_18519828/,"21"); print}' shunda_snps_filtered_temp14.vcf > shunda_snps_filtered_temp15.vcf
awk '{gsub(/PGA_scaffold8__67_contigs__length_18086311/,"10"); print}' shunda_snps_filtered_temp15.vcf > shunda_snps_filtered_temp16.vcf
awk '{gsub(/PGA_scaffold6__72_contigs__length_16467905/,"13"); print}' shunda_snps_filtered_temp16.vcf > shunda_snps_filtered_temp17.vcf
awk '{gsub(/PGA_scaffold5__81_contigs__length_16821318/,"14"); print}' shunda_snps_filtered_temp17.vcf > shunda_snps_filtered_temp18.vcf
awk '{gsub(/PGA_scaffold4__53_contigs__length_16101261/,"3"); print}' shunda_snps_filtered_temp18.vcf > shunda_snps_filtered_temp19.vcf
awk '{gsub(/PGA_scaffold3__61_contigs__length_14436489/,"5"); print}' shunda_snps_filtered_temp19.vcf > shunda_snps_filtered_temp20.vcf
awk '{gsub(/PGA_scaffold23__90_contigs__length_32925445/,"13"); print}' shunda_snps_filtered_temp20.vcf > shunda_snps_filtered_temp21.vcf
awk '{gsub(/PGA_scaffold1__43_contigs__length_10481396/,"1"); print}' shunda_snps_filtered_temp21.vcf > shunda_snps_filtered_temp22.vcf
awk '{gsub(/PGA_scaffold2__66_contigs__length_13522551/,"18"); print}' shunda_snps_filtered_temp22.vcf > shunda_snps_filtered_temp23.vcf


cp shunda_snps_filtered_temp23.vcf shunda_snps_filtered.vcf
rm shunda_snps_filtered_temp*



