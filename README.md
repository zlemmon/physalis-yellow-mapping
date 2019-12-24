# whole genome illumina alignment and variant calling

This repository contains code utilizing Docker and the Snakemake pipeline language to help build reproducible analyses. The expected inputs are a collection of Whole Genome Sequencing (WGS) Illumina reads and a reference genome. The output will be duplicates marked BAM files, unfiltered VCF files, and filtered VCF files by a fairly conservative cutoff (biallelic SNPs only excluding sites with indels with 100bp and with a map quality "MQ" >= 50). 

## Requisite software

A working installation of Docker and that's it

## Running the package

The Docker container is currently set up to run in interactive mode (`-it`). In future may convert to an auto-run. 

```code
docker run --rm -it -v ${PWD}/:/app bsa
```

Inside the container you will likely have to invoke `bwa index` on the reference fasta file to build the BWA indices (~15-20m for a 1-2 Gbp genome, but depends on how large your reference genome). Also make sure you point the `config*json` at the appropriate reference. When ready to run the main analysis do the following:

```code
snakemake -j --configfile /path/to/config.json
```
