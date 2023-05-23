#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
                     }

# replace file extension to "names"
outname <- gsub("(\\.[A-Za-z0-9]+)$", "\\.names", args[1])

res <- read.table(args[1], header=TRUE, row.names=NULL,sep = "\t")
res <- as.data.frame(res)

suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(dplyr))

# read genome annotation, it should coincide with RSEM/STAR annotation used
gr  <- read_gff("/usr/local/genomes/hg38.gtf")
gr  <- filter (gr, type == "gene" ) %>% select(gene_id, gene_name)
out <- as.data.frame(gr)

# only gene IDs and names
out <- out %>% select (gene_id, gene_name)
out <- out %>% mutate (gene_id = gsub("\\.\\d+$","",gene_id))

# replace NAs in gene_name to gene_ids
out <- out %>% mutate (gene_name = coalesce(gene_name, gene_id))

# merging tables
res.sub <- merge (res, out, by.x="GeneID", by.y="gene_id", all.x = TRUE, all.y = FALSE)

# replace NAs in gene_name of merged table to gene_ids
res.sub <- res.sub %>% mutate (gene_name = coalesce(gene_name, GeneID))

#                         new name    old name
res.sub = rename(res.sub, GeneName = gene_name)

# reorder columns, if move column to end than use like -GeneID
res.sub <- res.sub %>% select(GeneID, GeneName, everything())

write.table(res.sub, file=outname, quote = FALSE, sep='\t', row.names=FALSE)
