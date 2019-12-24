from ubuntu:16.04

ENV PATH=/root/miniconda3/bin:${PATH}

RUN apt-get update && apt-get install -y \
	wget

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
	mkdir /root/.conda && \
	bash Miniconda3-latest-Linux-x86_64.sh -b && \
	conda config --add channels defaults && \
	conda config --add channels bioconda && \
	conda config --add channels conda-forge

RUN conda install -y \
	bcftools \
	bwa \
	htslib \
	multiqc \
	picard \
	samtools \
	snakemake \
	trim-galore \
	;

