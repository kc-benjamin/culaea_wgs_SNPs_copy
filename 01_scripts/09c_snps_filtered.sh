#!/bin/bash
#SBATCH --partition=batch
#SBATCH --partition=batch
#SBATCH --job-name="09c_snps_filtered"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:20:00
#SBATCH --mem=1G
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=09b_plink_reformatting_%j.out


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
