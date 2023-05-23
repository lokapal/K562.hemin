#!/bin/sh
featureCounts -p -a /usr/local/genomes/hg38.gtf -o counts.txt -T 12 -t gene -g gene_id -M --readExtension5 2500 --readExtension3 2500 \
..\Align\rep1.white.bam ..\Align\rep2.white.bam ..\Align\rep1.red.bam ..\Align\rep2.red.bam
perl ./lib/counts_mselect.pl
