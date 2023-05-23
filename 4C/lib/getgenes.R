#!/usr/bin/Rscript
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(dplyr))
gr  <- read_gff("/usr/local/genomes/hg38.gtf")
gr  <- filter (gr, type == "gene" ) %>% select(gene_id, gene_name)
out <- as.data.frame(gr)
out <- out %>% select (seqnames, start, end, gene_name, gene_id, strand)
# replace NAs in gene_name to gene_ids
out <- out %>% mutate (gene_name = coalesce(gene_name, gene_id))
write.table(out, file="hg38.genes.bed", row.names=F, col.names=F, sep="\t", quote=F)
