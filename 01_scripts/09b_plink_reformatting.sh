#!/bin/bash
#SBATCH --partition=batch
#SBATCH --partition=batch
#SBATCH --job-name="09b_plink_reformatting"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:20:00
#SBATCH --mem=1G
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=09b_plink_reformatting_%j.out

sed 's/PGA_scaffold22__123_contigs__length_35631191/22/g' tmp2.bim > tmp3.bim
sed 's/PGA_scaffold21__64_contigs__length_27087477/21/g' tmp3.bim > tmp4.bim
sed 's/PGA_scaffold20__81_contigs__length_24420605/20/g' tmp4.bim > tmp5.bim
sed 's/PGA_scaffold18__88_contigs__length_23069817/18/g' tmp5.bim > tmp6.bim
sed 's/PGA_scaffold19__71_contigs__length_22530169/19/g' tmp6.bim > tmp7.bim
sed 's/PGA_scaffold17__83_contigs__length_22236360/17/g' tmp7.bim > tmp8.bim
sed 's/PGA_scaffold14__88_contigs__length_21401847/14/g' tmp8.bim > tmp9.bim
sed 's/PGA_scaffold16__56_contigs__length_21394868/16/g' tmp9.bim > tmp10.bim
sed 's/PGA_scaffold15__61_contigs__length_21329229/15/g' tmp10.bim > tmp11.bim
sed 's/PGA_scaffold12__56_contigs__length_20552211/12/g' tmp11.bim > tmp12.bim
sed 's/PGA_scaffold10__89_contigs__length_20417627/10/g' tmp12.bim > tmp13.bim
sed 's/PGA_scaffold13__42_contigs__length_19371112/13/g' tmp13.bim > tmp14.bim
sed 's/PGA_scaffold11__71_contigs__length_19041060/11/g' tmp14.bim > tmp15.bim
sed 's/PGA_scaffold7__70_contigs__length_18533579/7/g' tmp15.bim > tmp16.bim
sed 's/PGA_scaffold9__59_contigs__length_18519828/9/g' tmp16.bim > tmp17.bim
sed 's/PGA_scaffold8__67_contigs__length_18086311/8/g' tmp17.bim > tmp18.bim
sed 's/PGA_scaffold5__81_contigs__length_16821318/5/g' tmp18.bim > tmp19.bim
sed 's/PGA_scaffold6__72_contigs__length_16467905/6/g' tmp19.bim > tmp20.bim
sed 's/PGA_scaffold4__53_contigs__length_16101261/4/g' tmp20.bim > tmp21.bim
sed 's/PGA_scaffold3__61_contigs__length_14436489/3/g' tmp21.bim > tmp22.bim
sed 's/PGA_scaffold23__90_contigs__length_32925445/23/g' tmp22.bim > tmp23.bim
sed 's/PGA_scaffold2__66_contigs__length_13522551/2/g' tmp23.bim > final_bim_file.bim
sed 's/PGA_scaffold1__43_contigs__length_10481396/1/g' final_bim_file.bim > final_bim_file_fixed.bim

rm tmp*.bim final_bim_file.bim
mv final_bim_file_fixed.bim MU_snps_geno20_maf01_pruned_FIXED.bim
mv MU_snps_geno20_maf01_pruned.bim MU_snps_geno20_maf01_pruned_VER2.bim
mv MU_snps_geno20_maf01_pruned_VER2.bim 97_for_holding
mv MU_snps_geno20_maf01_pruned_FIXED.bim MU_snps_geno20_maf01_pruned.bim


