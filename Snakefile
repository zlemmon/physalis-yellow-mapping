"""
Author: Zachary H. Lemmon
Affiliation: Independent
Email: zlemmon@live.com
Aim: Snakemake workflow to process Whole Genome Sequencing (WGS) paired-end reads.
Date: August 26, 2019
Run: snakemake -s Snakefile
"""

# Globals ----------------------------------------------------
configfile: 'config.json'

FASTQFILES = config['fastqFiles']
SAMPLES = list(set([re.sub("_R[12]_001.fastq.gz", "", x) for x in FASTQFILES]))
REFERENCE = config['reference']

# bash safe mode
shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")

# Rules ------------------------------------------------------

rule final:
	input:
		expand("cleaned/{sample}_R1_001_val_1.fq.gz", sample = SAMPLES),
		expand("cleaned/{sample}_R2_001_val_2.fq.gz", sample = SAMPLES),
		expand("cleaned/{sample}_R1_001.fastq.gz_trimming_report.txt", sample=SAMPLES),
		expand("cleaned/{sample}_R2_001.fastq.gz_trimming_report.txt", sample=SAMPLES),
		expand("fastqc/{sample}_R1_001_val_1_fastqc.html", sample=SAMPLES),
		expand("fastqc/{sample}_R2_001_val_2_fastqc.html", sample = SAMPLES),
		expand("fastqc/{sample}_R1_001_val_1_fastqc.zip", sample=SAMPLES),
		expand("fastqc/{sample}_R2_001_val_2_fastqc.zip", sample=SAMPLES),
		expand("logs/{sample}.flagstat.txt", sample = SAMPLES),
		expand("bam/{sample}.dups.bam", sample = SAMPLES),
		expand("bam/{sample}.dups.bam.bai", sample = SAMPLES),
		expand("logs/{sample}.dups.Metrics", sample = SAMPLES),
		expand("vcf/{sample}.vcf.gz", sample = SAMPLES),
		expand("vcf/{sample}.filter.vcf.gz", sample = SAMPLES),
		expand("vcf/{sample}.filter.txt.gz", sample = SAMPLES),
		expand("logs/{sample}.mpileup.log", sample = SAMPLES),
		

rule trim_galore:
	input:
		r1 = "data/{sample}_R1_001.fastq.gz",
		r2 = "data/{sample}_R2_001.fastq.gz",
	output:
		r1 = "cleaned/{sample}_R1_001_val_1.fq.gz",
		r2 = "cleaned/{sample}_R2_001_val_2.fq.gz",
		log1 = "cleaned/{sample}_R1_001.fastq.gz_trimming_report.txt",
		log2 = "cleaned/{sample}_R2_001.fastq.gz_trimming_report.txt",
	shell:
		"""
		trim_galore --gzip --output_dir cleaned/ --paired {input.r1} {input.r2}
		"""
		
rule fastqc:
	input:
		r1 = "cleaned/{sample}_R1_001_val_1.fq.gz",
		r2 = "cleaned/{sample}_R2_001_val_2.fq.gz",
	output:
		html1 = "fastqc/{sample}_R1_001_val_1_fastqc.html",
		html2 = "fastqc/{sample}_R2_001_val_2_fastqc.html",
		zip1 = "fastqc/{sample}_R1_001_val_1_fastqc.zip",
		zip2 = "fastqc/{sample}_R2_001_val_2_fastqc.zip",
	shell:
		"""
		mkdir -p fastqc
		fastqc -o fastqc {input.r1} {input.r2}
		"""

rule bwa:
	input:
		r1 = "cleaned/{sample}_R1_001_val_1.fq.gz",
		r2 = "cleaned/{sample}_R2_001_val_2.fq.gz",
		ref = REFERENCE,
	output:
		sortedbam = temp("bam/{sample}.sorted.bam"),
		sortedbai = temp("bam/{sample}.sorted.bam.bai"),
		bam = temp("bam/{sample}.bam"),
		flagstat = "logs/{sample}.flagstat.txt",
		dups="bam/{sample}.dups.bam",
		dupsbai="bam/{sample}.dups.bam.bai",
		log="logs/{sample}.dups.Metrics",
	threads: 8
	shell:
		"""
		####
		# -t # number of threads to use specified maximum above by threads
		# -M # set as supplemental alignments
		# -R # record the read group as this
		####

		bwa mem -t {threads} -M \
		-R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}_lib1\\tLB:lib1\\tPL:illumina" \
		{input.ref} \
		{input.r1} \
		{input.r2} \
		| samtools view -b -T {input.ref} - > {output.bam}

		samtools sort -o {output.sortedbam} {output.bam}
		samtools index {output.sortedbam}
			
		samtools flagstat {output.sortedbam} > {output.flagstat}
		picard -Xmx12g MarkDuplicates VALIDATION_STRINGENCY=LENIENT \
			$(printf 'INPUT=%s ' {output.sortedbam}) \
			OUTPUT={output.dups} METRICS_FILE={output.log} \
			VERBOSITY=WARNING 

		samtools index {output.dups}
		"""


rule snpcalls:
	input:
		bam="bam/{sample}.dups.bam",
		ref=REFERENCE,
	output:
		snps="vcf/{sample}.vcf.gz",
		log="logs/{sample}.mpileup.log",
	shell:
		"""
		bcftools mpileup --fasta-ref {input.ref} --max-depth 1000000 --max-idepth 1000000 \
			--annotate FORMAT/DP,FORMAT/AD \
			--no-BAQ --min-BQ 0 \
			-O b {input.bam} | bcftools call -mv -O z -o {output.snps} 2> {output.log}
		bcftools reheader -s newsamplenames {output.snps} > temp.vcf.gz
		mv -f temp.vcf.gz {output.snps}
		tabix {output.snps}
		"""


rule filter_snps:
	input:
		snps="vcf/{sample}.vcf.gz",
	output:
		filtered_snps="vcf/{sample}.filter.vcf.gz",
		filtered_tsv="vcf/{sample}.filter.txt.gz",
	shell:
		"""
		bcftools view -m2 -M2 {input.snps} | bcftools filter \
			--SnpGap 100 -i ' TYPE="SNP" & MQ>=50 ' -Oz -o {output.filtered_snps}
		tabix {output.filtered_snps}

		bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%DP[\\t%SAMPLE %GT %AD]\n' {output.filtered_snps} | awk 'BEGIN{{FS="\\t"}}\
		        {{printf("%s\\t%s\\t%s\\t%s\\t%s",$1, $2, $3, $4, $5)}}\
			{{for(i=6;i<=NF;i++){{\
				string=$i; split(string, SGD, " ");\
				sample=SGD[1] ; gt=SGD[2] ; ad=SGD[3];\
				split(ad,ad2,",");\
				printf("\\t%s %s %s %s", sample, gt, ad2[1], ad2[2]);\
			}}\
			printf("\\n");\
		}}' | gzip > {output.filtered_tsv}
		"""
