#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name="read_depth_samtools"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=00-12:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out
#SBATCH --error=98_log_files/%x_%j.err

# Load needed modules
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh 
conda activate /home/kcb95328/conda/envs/culaea_pkgs

# Global variables
POPULATION="/home/kcb95328/Info-Shunda" #change this line as needed
name=$(cut -f1 $POPULATION/SRR_Acc_List_SL.txt | sed -n "${SLURM_ARRAY_TASK_ID}p") #change as needed
file="/scratch/kcb95328/ShundaLakeBrooks/06_bam_files/${name}.trimmed.fastq.gz.sorted.bam"

OUTPUT="/scratch/kcb95328/ShundaLakeBrooks/11_read_depth" #change as needed; should be from the submit directory
#name=$(cut -f1 $POPULATION/SRR_Acc_List_SL.txt | sed -n "${SLURM_ARRAY_TASK_ID}p") #change as needed


#Calculate coverage for a specific chromosome viewable at the terminal level, requires no secondary analysis but is per-file:
#samtools coverage -r PGA_scaffold14__88_contigs__length_21401847 -w 113 --plot-depth --ascii $DATAFOLDER/${name}.trimmed.fastq.gz.sorted.bam > ${OUTPUT}/${name}_chr20_coverage.txt

#Calculate coverage for a chromosome to be used with R:
samtools depth -a -H -J -r PGA_scaffold14__88_contigs__length_21401847:11,712,677-19,202,098 -f $file -o $OUTPUT/${name}_chr20_11.7to19.2_read-depth.txt
#-a: output all positions, including those with zero coverage
#-H: show column names at the beginning of the output
#-J: include reads with deletions in depth computation

conda deactivate
