# *ALLSITE* Variant Calling from Illumina Reads using GPU

## Overview
This repository contains a GPU-accelerated Snakemake workflow that calls both variant and non-variant sites from multi-sample Illumina reads using the GATK4 GPU implementation in NVIDIA's Clara Parabricks container. The pipeline is compatible with tools such as [`pixy`](https://pixy.readthedocs.io/en/latest/). Users supply a sample sheet and a configuration YAML specifying project metadata, reference genome, output directory, container image, and pipeline parameters such as window size, number of intervals, and per-interval memory. Snakemake and Singularity must be installed, and runs should reserve ample CPU cores, GPU resources, and memory (≤80% of the host) because some rules require more than 20 threads and high memory usage.

## Structure and Key Components
- **Snakefile** – Imports the configuration and sample sheet, pulls in rule modules, and defines the final `all` target for QC reports and filtered VCFs.
- **Environment definitions** – Conda YAMLs provide toolchains: `envs/envs.yaml` for general bioinformatics utilities, `envs/gatk4.yaml` for GATK, and `envs/blast.yaml` for low-complexity filtering via BLAST’s `dustmasker`.
- **Rules** (in `rules/`):
  - **0.qc.smk** – FastQC on trimmed reads, BAM stats, windowed coverage, Qualimap, bcftools and vcftools statistics, and a final MultiQC summary.
  - **1.preprocessing.smk** – Trims raw FASTQ files with `fastp`, producing paired/unpaired reads and QC reports.
  - **2.mapping.smk** – Copies and indexes the reference, then maps reads using Clara Parabricks’ GPU-accelerated `pbrun fq2bam` command, embedding read-group information and writing BAM files.
  - **3.variant_calling.smk** – Runs GPU-enabled `pbrun haplotypecaller` per sample, splits the reference into N intervals, builds a GenomicsDB, genotypes all sites in parallel, and merges the interval VCFs.
  - **4.vcf_filtering.smk** – Filters out low-quality and low-complexity regions, derives reference calls and high-quality SNPs, then merges them into an all-site VCF ready for downstream analyses.
- **scripts/** – Includes a helper script to split a FASTA reference into evenly sized BED intervals, facilitating parallel joint genotyping.
- **test_data/** – Provides example FASTQ files, a sample sheet, and a small reference for demonstration.

## Getting Started
1. Prepare a sample sheet and edit `configuration/config.yaml` with project details and resource paths.
2. Ensure Snakemake and Singularity are installed (ideally by conda), then run the command shown in the [Usage](#usage) section with appropriate `--cores`, `--resources`, and memory limits.

## Files to prepare:
 - A sample sheet - sample_sheet.csv: a comma delimited file with 3 columns (no column name):
   - `sample`, `path_to_read1`, `path_to_read2`
 - Modify configuration file - `configuration/config.yaml`:
   - `project`: a name for your project
   - `reference`:  path to the reference fasta file
   - `sample_sheet`: path to the sample sheet prepared above
   - `outdir`: path to the output directory
   - `clara-parabricks`: path to the [clara-parabricks](https://docs.nvidia.com/clara/parabricks/latest/index.html) image, e.g. "docker://nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1"
   - `w_size`: non-overlapping window size for reporting sequencing depths across the genome
   - `split_n`: number of independent jobs of `gatk` for parallelism
   - `memory_gb_per_interval`: memory in Gb for each independent `gatk` job

## Filtering parameters
The following configuration options adjust variant filtering in `rules/4.vcf_filtering.smk`:

- `gq_min`: minimum genotype quality to keep a call (default `20`).
- `qual_min`: minimum site QUAL score (default `30`).
- `dp_max`: maximum average depth across samples (default `50`).
- `f_missing_max`: maximum allowed fraction of missing genotypes (default `0.5`).

## Usage:
`snakemake --use-conda --use-singularity --singularity-args '--nv -B .:/dum' --cores [ncpu] --resources gpus=[ngpu] mem_gb=[mem]`

## Notes:
 - `clara-parabricks` now require 38Gb GPU memory for `fq2bam`. Therefore, `--low-memory` option is used in this step.
 - In snakemake command line (see below), `[ncpu]` should be larger than 20 as all resource usages are already hardcoded and some of rules uses more than 20 cpus.
 - In snakemake command line (see below), `[mem]` should not be more than 80% of the physical memory.
   

