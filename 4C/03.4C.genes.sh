#!/bin/bash
# script to find intersections between repicates, to filter DFAM, and to assign intersections to genes
# and their expression values 
#
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
#
# Input:  1. table_hg38.rep1.txt replicate 1 table (output of the script 02.mapping.sh)
#         2. table_hg38.rep2.txt replicate 2 table (output of the script 02.mapping.sh changed to process rep2)
#         3. lib/K562.white.mean.TPM gene expression values from RNAseq output (e.g. output from RNASeq directory) 
# Output: 1. K562.4C.white.nodfam.intersect.txt        K562 white 4C genome mappings intersections rep1 and rep2, DFAM removed table
#         2. K562.4C.white.nodfam.intersect.bedGraph   K562 white 4C hg38 profile for genome browsers
#         3. K562.4C.white.intersect.all.txt           K562 white 4C intersect noDFAM with genes and expression added, all entries (with genes intersections and empty)
#         4. K562.4C.white.intersect.genes.txt         K562 white 4C intersect noDFAM with genes and expression added, genes intersections only
#
# Dependency tools:
# 1. bedtools 2.29.1      https://bedtools.readthedocs.io/en/latest/
# 2. BEDOPS 2.4.40        https://bedops.readthedocs.io/en/latest/
# 3. R with libraries refGenome and dplyr
# 4. DFAM database hg38   https://www.dfam.org/releases/Dfam_3.6/annotations/hg38/hg38_dfam.nrph.hits.gz
# 5. hg38 genome annotation in GTF format, unpacked, classic chromosome names like "chr1, chr2, ..."
#                         http://ftp.ensembl.org/pub/current_gtf/homo_sapiens/
#
# Part 1. Find intersection between replicates, remove mappings that are completely inside DFAM entries
# get DFAM database and convert it to bed format
# set type
TYPE=white
wget https://www.dfam.org/releases/current/annotations/hg38/hg38.nrph.hits.gz
./lib/dfam2bed.pl hg38.nrph.hits.gz
mv hg38.nrph.hits.bed hg38_dfam.bed
sort -k1,1 -k2,2n -k3,3n hg38_dfam.bed -o hg38_dfam.bed
#
cp table_hg38.rep1.txt rep1.txt
cp table_hg38.rep2.txt rep2.txt
./lib/table2bedg.pl rep1.txt
./lib/table2bedg.pl rep2.txt
sort -k1,1 -k2,2n -k3,3n rep1.bedGraph -o rep1.bedGraph
sort -k1,1 -k2,2n -k3,3n rep2.bedGraph -o rep2.bedGraph
bedtools intersect -sorted -v -f 1.0 -a rep1.bedGraph -b hg38_dfam.bed  > 4C_hg38_nodfam.rep1.bedGraph
bedtools intersect -sorted -v -f 1.0 -a rep2.bedGraph -b hg38_dfam.bed  > 4C_hg38_nodfam.rep2.bedGraph
bedtools intersect -a 4C_hg38_nodfam.rep1.bedGraph -b 4C_hg38_nodfam.rep2.bedGraph -wb > intersect_reps.txt
./lib/overlist2bed_mean.pl intersect_reps.txt
sort -k1,1 -k2,2n -k3,3n intersect_reps.bedGraph -o intersect_reps.bedGraph
sort -k1,1 -k2,2n -k3,3n intersect_reps.bed -o intersect_reps.bed
bedtools merge -c 4 -o mean -i intersect_reps.bed > intersect_merged.bed
./lib/mergebedgraphmean.sh intersect_reps.bedGraph K562.4C.$TYPE.nodfam.intersect.bedGraph
./lib/addSubseq.pl /usr/local/genomes/hg38.mfa intersect_merged.bed > K562.4C.$TYPE.nodfam.intersect.txt
sort -k1,1 -k2,2n -k3,3n K562.4C.$TYPE.nodfam.intersect.txt -o K562.4C.$TYPE.nodfam.intersect.txt
#
# Part 2. Assign genes to K562 4C-associated intersected replicates noDFAM mappings
cp K562.4C.$TYPE.nodfam.intersect.txt 4C.bed
# Annotation h38.gtf should be in /usr/local/genomes
Rscript ./lib/getgenes.R
sort -k1,1 -k2,2n -k3,3n hg38.genes.bed -o hg38.genes.bed
bedtools intersect -wa -wb -a 4C.bed -b hg38.genes.bed > 4C.genetable
bedtools intersect -v -wa -a 4C.bed -b hg38.genes.bed > 4C.empty
./lib/tbl_addexpr.pl 4C.genetable ./lib/K562.$TYPE.mean.TPM
cat 4C.genetable.expr 4C.empty > K562.4C.$TYPE.intersect.all.txt
mv 4C.genetable.expr K562.4C.$TYPE.intersect.genes.txt
sort -k1,1V -k2,2n -k3,3n K562.4C.$TYPE.intersect.genes.txt -o K562.4C.$TYPE.intersect.genes.txt
sort -k1,1V -k2,2n -k3,3n K562.4C.$TYPE.intersect.all.txt -o K562.4C.$TYPE.intersect.all.txt
