# *ALLSITE* Variant Calling from Illumina Reads using GPU

## Description:
 - A GPU-accelarated snakemake workflow that calls variants from multi-sample illumina reads using GATK4 GPU version.
 - It calls both variant and non-variant sites. Compatible with tools such as [`pixy`](https://pixy.readthedocs.io/en/latest/).

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

## Environment:
 - Make sure snakemake and singularity is installed in current environment.

## Notes:
 - `clara-parabricks` now require 38Gb memory for `fq2bam`. Therefore, `--low-memory` option is used.
 - In snakemake command line (see below), `[ncpu]` should be larger than 20 as all resource usages are already hardcoded and some of rules uses more than 20 cpus.
 - In snakemake command line (see below), `[mem]` should not be more than 80% of the physical memory.
   
## Usage:
`snakemake --use-conda --use-singularity --singularity-args '--nv -B .:/dum' --cores [ncpu] --resources gpus=[ngpu] mem_gb=[mem]`
