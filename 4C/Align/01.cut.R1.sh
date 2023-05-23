#!/bin/bash
TYPE=white
REP=rep1
# leading full adapters at 5' end cut
cutadapt -O 10 -j 20 --trim-n --times=4 --minimum-length 20 -q 24 --untrimmed-output=untr.R1.1.fastq.gz \
-g file:./lib/fulladapters.fa \
-o R1.1.fastq.gz 4C.K562.$TYPE.$REP.R1.fastq.gz
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
cat untr.R5.1.fastq.gz R5.1.fastq.gz  > R1.fastq.gz
#
# second step - we will try to find 4C adapters in UNTRIMMED reads
# now we require the more nucleotides from adapter (15), but higher alignment error rate (0.2)
cutadapt -O 15 -e 0.2 -j 20 --trim-n --times=4 --minimum-length 20 -q 24 --untrimmed-output=untr.untr.R1.1.fastq.gz \
-g file:./lib/fulladapters.all.fa \
-o full.LQ.R1.fastq.gz untr.R1.1.fastq.gz
#
cutadapt -O 5 -j 20 --times=4 --minimum-length 20 -q 24 --untrimmed-output=untr.R2.1.fastq.gz \
-a file:./lib/illumina3.adapters.fa \
-o R2.1.fastq.gz full.LQ.R1.fastq.gz
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
cat untr.R5.1.fastq.gz R5.1.fastq.gz  > R1.LQ.fastq.gz
cat R1.fastq.gz R1.LQ.fastq.gz > $TYPE.$REP.R1.fastq.gz
# cleanup from temporary files
rm -f untr*.fastq.gz R*.fastq.gz
