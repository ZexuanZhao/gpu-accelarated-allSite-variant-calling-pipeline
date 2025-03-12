# ALLSITE Variant Calling from Illumina Reads using GPU

## Description:
 - A GPU-accelarated snakemake workflow that calls variants from multi-sample illumina reads using Deepvariant and GLnexus

## Files to prepare:
 - A sample sheet - sample_sheet.csv: a comma delimited file with 3 columns (no column name):
   - `sample`, `path_to_read1`, `path_to_read2`
 - Modify configuration file - `configuration/config.yaml`:
   - `project`: a name for project
   - `reference`:  path to reference fasta file
   - `sample_sheet`: path to the sample sheet prepared above
   - `outdir`: path to the output directory
   - `clara-parabricks`: image path to clara-parabricks, e.g. "docker://nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1"
   - `w_size`: non-overlapping window size of reporting average depth along the genome.
   - `genome_size`: ref genome size in bp
   - `split_n`: number of independent jobs of `gatk` for parallelism
   - `memory_gb_per_interval`: memory in Gb for each independent `gatk` job
 
## Environment:
 - Make sure snakemake and singularity is installed in current environment.

## Notes:
 - `clara-parabricks` now require 38Gb memory for `fq2bam`. Therefore, `--low-memory` option is used.
 - In snakemake commandline, [ncpu] should be larger than 20 as all resource usages are hardcoded.
## Usage:
`snakemake --use-conda --use-singularity --singularity-args '--nv -B .:/dum' --cores [ncpu] --resources gpus=[ngpu]`
