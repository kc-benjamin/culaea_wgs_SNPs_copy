#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name="05_pat_on_back"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1K
#SBATCH --time=00-00:01:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out

echo "Give yourself a pat on the back. You've made it half way, and you're doing great!"
