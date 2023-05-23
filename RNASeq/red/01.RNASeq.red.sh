#!/bin/sh
# script to get and process gene expression in TPM of H.sapiens K562 incubated with 50 uM hemin for 108h hg38 
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Input:  H.sapiens K562 cells incubated with 50 uM hemin for 108h (red) 
#         RNA-Seq data: GEO GSE232390: GSM7329826, GSM7329827
# Output: 1. K562.red.RNAseq1.genes.results    K562 red hg38 gene expression replicate 1
#            K562.red.RNAseq2.genes.results    K562 red hg38 gene expression replicate 2
#         2. K562.red.hg38.TPM                 joined replicates expression all and average values
#            K562.red.hg38.names               joined replicates expression all and average values with ISO names
#            K562.red.hg38.mean.TPM            K562 red hg38 gene expression average values per gene
#            K562.red.hg38.mean.names          K562 red hg38 gene expression average values per gene with ISO names
#
# Dependency tools:
# 1. NCBI SRA Toolkit     https://github.com/ncbi/sra-tools
# 2. Trimmomatic 0.39     http://www.usadellab.org/cms/?page=trimmomatic
# 3. RSEM 1.3.1           https://github.com/deweylab/RSEM
# 4. STAR 2.7.5c          https://github.com/alexdobin/STAR
# 5. R with libraries     tibble, dplyr
# Get data from NCBI SRA Archive. XXXXXXX will be changed when data will be released to public
#fastq-dump --split-files --gzip SRRXXXXXXXX
#fastq-dump --split-files --gzip SRRXXXXXXXX
mv SRRXXXXX_1.fastq.gz red.rep1.fastq.gz
mv SRRXXXXX_1.fastq.gz red.rep2.fastq.gz
# Quality trimming
trimmomatic SE -threads 30 red.rep1.fastq.gz rep1.filtered.fastq.gz LEADING:18 TRAILING:18 SLIDINGWINDOW:4:22 MINLEN:20 ILLUMINACLIP:TruSeq3-SE:2:30:10
trimmomatic SE -threads 30 red.rep2.fastq.gz rep2.filtered.fastq.gz LEADING:18 TRAILING:18 SLIDINGWINDOW:4:22 MINLEN:20 ILLUMINACLIP:TruSeq3-SE:2:30:10
# Expression calculation per replicate
rsem-calculate-expression -p 32 --fragment-length-mean 255 --star -p 32 --output-genome-bam --calc-ci \
--ci-memory 30720 --star-gzipped-read-file rep1.filtered.fastq.gz /usr/local/genomes/hg38.rsem/hg38 K562.red.RNAseq1
rsem-calculate-expression -p 32 --fragment-length-mean 255 --star -p 32 --output-genome-bam --calc-ci \
--ci-memory 30720 --star-gzipped-read-file rep2.filtered.fastq.gz /usr/local/genomes/hg38.rsem/hg38 K562.red.RNAseq2
# Calculate average expression values between replicas
./lib/mean_expr.R
./lib/addisonames.R K562.red.hg38.TPM
./lib/addisonames.R K562.red.hg38.mean.TPM
rm -f *.bam
