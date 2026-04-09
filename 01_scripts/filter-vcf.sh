#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name="filter-vcf"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00-12:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out
#SBATCH --error=98_log_files/%x_%j.err

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh 
conda activate /home/kcb95328/conda/envs/culaea_pkgs

#bcftools index muir_snps_filtered.vcf.gz

bcftools view --regions 20 shunda_snps_filtered.vcf.gz > shunda_snps_filtered_chr20.vcf

conda deactivate