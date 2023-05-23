#!/bin/sh
# Dependencies: featureCounts      https://subread.sourceforge.net/
featureCounts -a /usr/local/genomes/hg38.gtf -o counts.txt -T 30 -t exon -g gene_id \
../white/rnaseq.white.rep1.bam ../white/rnaseq.white.rep2.bam \
../red/rnaseq.red.rep1.bam ../red/rnaseq.red.rep2.bam
