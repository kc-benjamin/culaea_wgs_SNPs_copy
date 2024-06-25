# culaea_wgs_SNPs
This directory includes the scripts used to call SNPs from raw fastq files generated from Illumina whole-genome shotgun sequence data.

The scripts are in a subdirectory called 01_scripts, and are numbered 01-09. They were written to run on a SLURM-based HPC, and some scripts have companion BASH scripts to run parallel jobs for each sample (01, 02, 03, 04, and 06) or for each chromosome (script 07).

Scripts 05, 08, and 09 are run as single jobs.

More details on these scripts are below.

The scripts expect a directory structure that can be created with the following command.
```
mkdir  01_scripts  02_info_files  03_genome  04_raw_data  05_trimmed_data  06_bam_files  07_raw_VCFs  98_log_files  99_metrics
```

The 02_info_files directory should contain the following files:
* `datatable.txt` - contains all the metadata for each individual in the dataset, with the sample name for each individual in the first column. An example is provided.
* `sample_file.txt` - just a list (in a single column) of the sample names for all individuals (which might seem redundant, and it is). An example is provided.
* `chromosomes.txt` - a list (in a snigle column) of all the chromosomes in the reference genome (can be created form the .dict file that should be in the 03_genome directory). An example is provided.
* `bamfiles.txt` - a list (in a single column) of the paths to all the bam files generated in step 06, that will be used in step 07 (i.e. all the "realigned.bam" files). An example is provided.

The 03_genome directory should contain an indexed reference genome. A script to index the reference genome is included here. IMPORTANT NOTE: for some reason, the extension on your reference genome has to be .fa (not .fna or anything else) for this set of scripts to work!!!

The 04_raw_data directory should contain all the raw fastq files (or fastq.gz files) that you want to analyze.

All the other directories should start off empty.

## Step 1 - Read trimming
Uses fastp.

* Scripts: `01_fastp.sh` and `01_fastp_submit`
* Inputs: raw fastq (or fastq.gz) files
* Outputs: trimmed fastq files

## Step 2 - Alignment to reference genome
Uses bwa-mem. Also adds read group info.

* Scripts: `02_bwa_alignments.sh` and `02_bwa_alignments_submit`
* Inputs: trimmed fastq files
* Outputs: sorted and indexed bam files with read group information in header lines

## Step 3 - Get BAM quality metrics
Uses Picard

* Scripts: `03_collect_metrics.sh` and `03_collect_metrics_submit`
* Inputs: sorted and indexed bam files
* Outputs: alignment and insert size metrics

## Step 4 - Remove duplicates
Uses Picard

* Scripts: `04_remove_duplicates.sh` and `04_remove_duplicates_submit`
* Inputs: sorted and indexed bam files
* Outputs: Deduplicated BAMs

## Step 6 - Realignment around indels
Uses GATK's `RealignerTargetCreator`

* Scripts: `06_gatk_realignments.sh` and `06_gatk_realignments_submit`
* Inputs: Deduplicated BAMs
* Outputs: Realigned BAMs

## Step 7 - SNP calling
Uses mpileup (in bcftools)

* Scripts: `07_bcftools_vc.sh` and `07_bcftools_vc_submit`
* Inputs: Realigned BAMs (all of them - make sure step 6 is complete for all samples before starting this step)
* Outputs: Unfiltered VCFs

## Step 8 - Concatenate and filter
Uses bcftools and vcftools

* Script: `08_concat.sh`
* Inputs: Unfiltered VCFs
* Outputs: a filtered VCF file (`*_filtered.vcf.gz`, where * denotes a descriptive name for your dataset that you specify in this script)

## Step 9 - further SNP filtering and file processing
Uses PLINK and vcftools and GCTA. This script is designed to accommodate non-model organisms. You'll have to modify it to reflect the names of the chromosomes (or chromosome-level contigs) in your reference genome. This is because PLINK (and other programs) only know how interpret particular chromosome names (i.e. 1, 2, 3, etc.).

* Script: `09_process_snps.sh`
* Inputs: `*_filtered.vcf.gz`
* Outputs: plink-formatted version of your SNPs filtered for missingness and minimum allele frequency, LD-pruned SNPs, and a GRM.

