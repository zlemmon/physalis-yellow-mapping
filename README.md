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

This gets you three separate VCF files and alignment bam files. I further processed by filtering the Phy50 yellow parent to only include homozygous sites and then recalling SNPs for the genome across 3 alignments all at once. 

```bash
bcftools filter -i 'GT=="1/1"' vcf/Phy50_Lib_S3.filter.vcf.gz | bcftools view -Oz -o vcf/Phy50_Lib_S3.filter.homozygous.vcf.gz

bcftools mpileup --fasta-ref ref/G20_hifi_hybrid_scaffold_ncbi_notScaffolded.fasta --max-depth 1000000 --max-idepth 1000000 --annotate FORMAT/DP,FORMAT/AD --no-BAQ --min-BQ 0 -T vcf/Phy50_Lib_S3.filter.homozygous.vcf.gz -O b bam/Phy50_Lib_S3.dups.bam bam/purple_pool_S1.dups.bam bam/yellow_pool_S2.dups.bam | bcftools call -mv -O z -o vcf/combined_calls.vcf.gz
```

At this point, we should have the combined variant calls across all three libraries in the `vcf/combined_calls.vcf.gz` file. They can be parsed into a gunzip compressed text tab-delimited file with the following rerun of the Snakefile

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP[\t%SAMPLE %GT %AD]\n' vcf/combined_calls.vcf.gz | awk 'BEGIN{FS="\t"} \
	{printf("%s\t%s\t%s\t%s\t%s",$1, $2, $3, $4, $5)} \
	{for(i=6;i<=NF;i++){ \
		string=$i; split(string, SGD, " "); \
		sample=SGD[1] ; gt = SGD[2] ; ad = SGD[3]; \
		split(ad, ad2, ","); \
		printf("\t%s %s %s %s", sample, gt, ad2[1], ad2[2]); \
	} \
	printf("\n"); }' | gzip > vcf/combined_calls.txt.gz
```
