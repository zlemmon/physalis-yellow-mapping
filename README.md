# whole genome illumina alignment and variant calling

This repository contains code utilizing Docker and the Snakemake pipeline language to help build reproducible analyses. The expected inputs are a collection of Whole Genome Sequencing (WGS) Illumina reads (in the `data/` directory) and a reference genome fasta file (in the `ref/` directory). The output will be duplicates marked BAM files, unfiltered VCF files, and filtered VCF files by a fairly conservative cutoff (biallelic SNPs only excluding sites with indels with 100bp and with a map quality "MQ" >= 50). 

## Requisite software

A working installation of Docker and that's it

## Running the package

A starting `config.json` file is needed that has two key-value pairs, `fastqFiles` should be a json list of fastq files and `reference` should be a string with the full path to the reference fasta genome file. Place your raw read fastq files in the `data/` directory. You will need to add the read file names to the `fastqFiles` list in `config.json`. Please add both read1 and read2 (do not include the data/ prefix). You will also need to place your reference genome fasta file in the `ref/` directory and add this file to the `reference` key in `config.json` (include the `ref/` prefix). All other analysis directories should be automatically created during the Snakemake pipeline run as needed. 

The Docker container is currently set up to run in interactive mode (`-it`). In future may convert to an auto-run Docker entrypoint, but for now since we assume the bwa indices may not be pre-generated we run interactively. 

```code
docker run --rm -it -v ${PWD}/:/app bsa
```

Inside the container you will likely have to invoke `bwa index` on the reference fasta file to build the BWA indices (~15-30m for a 1-2 Gbp genome, but depends on how large your reference genome is and how big a machine you can throw at it). Also make sure you point the `config*json` at the appropriate reference. When ready to run the main analysis do the following from in the running container:

```code
snakemake -j --configfile /path/to/config.json
```
