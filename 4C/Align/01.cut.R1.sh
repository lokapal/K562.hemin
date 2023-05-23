#!/bin/bash
TYPE=white
REP=rep1
# leading full adapters at 5' end cut
cutadapt -O 10 -j 20 --trim-n --times=4 --minimum-length 20 -q 24 --untrimmed-output=untr.R1.1.fastq.gz \
-g file:./lib/fulladapters.fa \
-o R1.1.fastq.gz 4C.K562.minus.rep1.R1.fastq.gz
mv R1.1.fastq.gz R1.j1.fastq.gz
#cat untr.R1.1.fastq.gz 
#
cutadapt -O 5 -j 20 --times=4 --minimum-length 20 -q 24 --untrimmed-output=untr.R2.1.fastq.gz \
-a file:./lib/illumina3.adapters.fa \
-o R2.1.fastq.gz R1.j1.fastq.gz
cat untr.R2.1.fastq.gz R2.1.fastq.gz  > R2.j1.fastq.gz
#
cutadapt -O 10 -j 20 --times=2 --minimum-length 20 -q 24 --untrimmed-output=untr.R3.1.fastq.gz \
-g file:./lib/fulladapters5.all.fa \
-o R3.1.fastq.gz R2.j1.fastq.gz
cat untr.R3.1.fastq.gz R3.1.fastq.gz  > R3.j1.fastq.gz
#
cutadapt -O 5 -j 20 --times=4 --minimum-length 20 -q 24 --untrimmed-output=untr.R4.1.fastq.gz \
-a file:./lib/illumina.adapters.fa \
-o R4.1.fastq.gz R3.j1.fastq.gz
cat untr.R4.1.fastq.gz R4.1.fastq.gz  > R4.j1.fastq.gz
#
cutadapt -O 10 -j 20 --times=4 --minimum-length 20 -q 24 --untrimmed-output=untr.R5.1.fastq.gz \
-a file:./lib/fulladapters.tails.fa \
-o R5.1.fastq.gz R4.j1.fastq.gz
cat untr.R5.1.fastq.gz R5.1.fastq.gz  > R5.j1.fastq.gz
