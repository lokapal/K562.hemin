#!/bin/sh
# script to align RNASeq data with the same parameters as in rsem 1.3.1
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Input:  pre-filtered RNASeq rep1 and RNASeq rep2 datasets
# Output: 1. rnaseq.white.rep1.bam          K562 untreated RNASeq replicate 1 alignment to hg38
#            rnaseq.white.rep2.bam          K562 untreated RNASeq replicate 2 alignment to hg38
#         2. WHITE_Pearson.pdf              K562 untreated (white) RNASeq replicates consistency check: Pearson correlation
#            WHITE_Spearman.pdf             K562 untreated (white) RNASeq replicates consistency check: Spearman correlation
#
# Dependency tools:
# 1. STAR 2.7.5c          https://github.com/alexdobin/STAR
# 2. samtools 1.15        http://www.htslib.org
# 3. deepTools2 3.5.1     https://deeptools.readthedocs.io
# align to hg38 genome with the options that are the same as in rsem 1.3.1
STAR --genomeDir /usr/local/genomes/hg38.rsem --outSAMunmapped Within --outFilterType BySJout \
--outSAMattributes NH HI AS NM MD --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --runThreadN 32 \
--genomeLoad NoSharedMemory --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM \
--outFileNamePrefix rep1. --readFilesCommand zcat --readFilesIn rep1.filtered.fastq.gz
# 
STAR --genomeDir /usr/local/genomes/hg38.rsem --outSAMunmapped Within --outFilterType BySJout \
--outSAMattributes NH HI AS NM MD --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --runThreadN 32 \
--genomeLoad NoSharedMemory --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM \
--outFileNamePrefix rep2. --readFilesCommand zcat --readFilesIn rep2.filtered.fastq.gz
mv rep1.Aligned.sortedByCoord.out.bam rnaseq.white.rep1.bam
samtools index -@ 16 rnaseq.white.rep1.bam
mv rep2.Aligned.sortedByCoord.out.bam rnaseq.white.rep2.bam
samtools index -@ 16 rnaseq.white.rep2.bam
# Quality Control to check consistency between replicates (optional)
bamCoverage --effectiveGenomeSize 2913022398 --normalizeUsing RPKM --exactScaling -p 30 -b rnaseq.white.rep1.bam -o rep1.bw
bamCoverage --effectiveGenomeSize 2913022398 --normalizeUsing RPKM --exactScaling -p 30 -b rnaseq.white.rep2.bam -o rep2.bw
multiBigwigSummary bins -p 30 -b rep1.bw rep2.bw -o res12hs.npz
plotCorrelation -in res12hs.npz -c spearman --removeOutliers --skipZeros -p scatterplot -o WHITE_spearman.pdf --labels WHITE.rep1 WHITE.rep2 --log1p
plotCorrelation -in res12hs.npz -c pearson  --removeOutliers --skipZeros -p scatterplot -o WHITE_pearson.pdf  --labels WHITE.rep1 WHITE.rep2 --log1p
