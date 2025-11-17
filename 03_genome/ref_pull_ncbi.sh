#!/bin/bash
#SBATCH --partition=batch
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=ref_pull_ncbi_%j.out
#SBATCH --job-name=ref_index
#SBATCH --error=ref_pull_ncbi_%j.err

eval "$(conda shell.bash hook)"
conda activate /home/kcb95328/miniconda3/envs/ncbi_datasets
conda activate ncbi_datasets
    conda install -c conda-forge ncbi-datasets-cli

module load datasets/3.1.0-foss-2023a

datasets download genome accession GCF_949316345.1 --filename punpun-ref-genome.zip
unzip punpun-ref-genome.zip -d 03_genome/