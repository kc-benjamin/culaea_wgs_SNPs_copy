#!/bin/bash

#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=4G
#SBATCH --time=0-12:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out
#SBATCH --error=98_log_files/%x_%j.err
#SBATCH --array=1-23

# Load needed modules
#module load BCFtools/1.21-GCC-13.3.0
#gcc/9.3.0
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh 
conda activate /home/kcb95328/conda/envs/culaea_pkgs

# Global variables
POPULATION="/home/kcb95328/Info-Muir"
GENOMEFOLDER="/home/kcb95328/genomes"
GENOME=$(ls -1 $GENOMEFOLDER/brook_genome_hap1_v1.fasta | xargs -n 1 basename)
VCF="07_vcfs"
BAM="ML_bamfiles.txt" #this should be generated at the command line from the work dir each time this is run: find . -name "*trimmed.fastq.gz.sorted.bam" -printf "%h/%f\n" > ML_bamfiles.txt
echo $BAM
SAMPS="$POPULATION/SRR_Acc_List_ML.txt"

#Pass the chromosome number from sbatch command
chrom_num=$SLURM_ARRAY_TASK_ID

# Fetch chromosome from the array
CHROM=$(sed -n "${chrom_num}p" $GENOMEFOLDER/brook_genome_hap1_v1_chromosomes2.txt)
echo $CHROM
#CHROM="PGA_scaffold2__66_contigs__length_13522551"
#SCAFFOLD=$(echo "$CHROM" | grep -oP 'scaffold\d+')

#removed the -q 5 quality filtering from this step
bcftools mpileup -Ou --fasta-ref $GENOMEFOLDER/$GENOME --bam-list $BAM -r $CHROM -I -a FMT/AD | \
	bcftools call -S "$SAMPS" -G - -f GQ -mv -Ov > "$VCF/${CHROM}.vcf"

conda deactivate