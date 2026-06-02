#!/bin/bash
#SBATCH --partition=batch
#SBATCH --partition=batch
#SBATCH --job-name="featureCounts"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00
#SBATCH --mem=1G
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=featureCounts_%j.out
#SBATCH --error=featureCounts_%j.err
#SBATCH --array=1-40

ml Subread/2.0.6-GCC-12.3.0

#variables
INDIR="03_star_output"
OUTDIR="04_counts"
GENOMEDIR="02_genome"

#Pass the sample number from the sbatch command
samp_num=$SLURM_ARRAY_TASK_ID

# Pull sample name from the sample info
name=$(cut -f1 raw_data_list3.txt | sed -n "${samp_num}p")

subread featureCounts -T $SLURM_CPUS_PER_TASK \
     -a $GENOMEDIR/brook_genome_hap1_v1.gff \
     --input_files $INDIR/${name}Aligned.sortedByCoord.out.bam \
     -o $OUTDIR/${name}_counts.txt 