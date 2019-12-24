#!/bin/bash

# trim reads
trim_galore --gzip --output_dir cleaned/ --paired data/Phy50_Lib_S3_R1_001.fastq.gz data/Phy50_Lib_S3_R2_001.fastq.gz

#*file

#bwa mem -1 $file1 -2 $file2 | samtools view -b > out.vcf.gz

#bcftools view stuff ....
