#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name="read_depth"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=00-12:00:00
#SBATCH --mail-user=kcb95328@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=98_log_files/%x_%j.out
#SBATCH --error=98_log_files/%x_%j.err

#load necessary modules
ml deepTools/3.5.5-gfbf-2023a

#variables:
POPULATION="/home/kcb95328/Info-Shunda"
#samp_num=$SLURM_ARRAY_TASK_ID
#name=$(cut -f1 $POPULATION/SRR_Acc_List_SL.txt | sed -n "${samp_num}p") #change as needed

 # 1. Loop through all BAMs and create bigwigs for chr1
for bam in 01_aligned_bams/*.trimmed.fastq.gz.sorted.bam; do
    echo "Processing $bam..."
    bamCoverage -b "$bam" -o "${bam%.trimmed.fastq.gz.sorted.bam}_chr1.bw" --region PGA_scaffold14__88_contigs__length_21401847 --binSize 100 --normalizeUsing RPKM --numberOfProcessors 8
done

#2: Compile the data matric with deeptools
multiBigwigSummary bins \
    -b *_chr20.bw \
    -o all-bams_chr20.npz \
    --region PGA_scaffold14__88_contigs__length_21401847 \
    --binSize 100

#3: Generate the scan
plotProfile \
    -m all-bams_chr20.npz \
    -o all-bams_chr20_profile.png \
    --plotTitle "Normalized Read Depth Across Chromosome 20 for all samples" \
    --regionsLabel "Chromosome 20" \
    --plotType lines
