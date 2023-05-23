#/bin/bash
TYPE=white
REP=rep1
repair.sh -Xmx50g in1=$TYPE.$REP.R1.fastq.gz in2=$TYPE.$REP.R2.fastq.gz \
out1=4C.$TYPE.$REP.filtered.R1.fastq.gz out2=4C.$TYPE.$REP.filtered.R2.fastq.gz \
outs=4C.K562.filtered.singletons.fastq.gz repair overwrite=t
