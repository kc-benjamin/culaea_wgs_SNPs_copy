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

bcftools view --regions PGA_scaffold14__88_contigs__length_21401847:15,437,889-15,458,061 muir_snps_filtered.vcf.gz > muir_snps_filtered_amhy.vcf

conda deactivate