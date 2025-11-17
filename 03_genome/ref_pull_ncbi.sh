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

conda init bash
conda activate /home/kcb95328/miniconda3/envs/ncbi_datasets
conda activate ncbi_datasets
    conda install -c conda-forge ncbi-datasets-cli

module load datasets/3.1.0-foss-2023a

datasets download genome accession CASGFK000000000.1 --filename punpun-ref-genome.fa