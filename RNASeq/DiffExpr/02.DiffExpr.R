#!/usr/bin/Rscript
# script to perform differential RNASeq analysis 
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Input:  counts.txt pre-computed file by featureCounts
# Output: 1. RNASEq.WHITE-RED.results.tsv          K562 White-Red Differential expression analysis tabbed file
#         2. RNASEq.WHITE-RED.scatterplots.pdf     Scatterplots to control replicates consistency
#         3. RNASEq.WHITE-RED.volcanoplot.pdf      VolcanoPlot to display the most prominent down/upregulated genes
#         4. RNASEq.WHITE-RED.PCA.pdf              PCA to control replicates consistency
#
# Dependency tools & libraries:
# 1. R
# Bioconductor packages:
# 2. DESeq2 R             https://bioconductor.org/packages/release/bioc/html/DESeq2.html
# 3. EnhancedVolcano      https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html
# 4. rtracklayer          https://bioconductor.org/packages/release/bioc/html/rtracklayer.html
# 5. genefilter           https://bioconductor.org/packages/release/bioc/html/genefilter.html
# 6. PCAExplorer          https://bioconductor.org/packages/release/bioc/html/pcaExplorer.html
# Common R libraries:
# 6. dplyr, ggplot2, tibble, RColorBrewer, gplots, ggrepel, calibrate

# Import data from featureCounts
## Previously ran at command line something like this:
## featureCounts -a genes.gtf -o counts.txt -T 12 -t exon -g gene_id GSM*.sam
countdata <- read.table("counts.txt", header=TRUE, row.names=1)

# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,6:ncol(countdata)]

# Remove .bam or .sam from filenames
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))
colnames(countdata) <- gsub("\\.\\.\\.", "", colnames(countdata))

# Convert to matrix
countdata <- as.matrix(countdata)

# Assign condition (first pair is control, second pair is experiment)
(condition <- factor(c(rep("ctl", 2), rep("exp", 2))))

# Analysis with DESeq2
suppressPackageStartupMessages(library(DESeq2))

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)

vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
rld <- rlogTransformation(dds)
hist(assay(rld))

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tibble))

ddc <- estimateSizeFactors(dds)

df <- bind_rows(
  as_tibble(log2(counts(ddc, normalized=TRUE)[, 1:2]+1)) %>% mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "VarianceStabilizingTransformation"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rLogTransformation (DESeq2 default)"))

colnames(df)[1:2] <- c("x", "y")  

cairo_pdf("RNASEq.WHITE-RED.scatterplots.pdf",width=15,height=10,antialias="default",fallback_resolution = 300,onefile=T)
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation) + labs(x = "rep1.WHITE", y="rep2.WHITE")

rm (df)
df <- bind_rows(
  as_tibble(log2(counts(ddc, normalized=TRUE)[, 3:4]+1)) %>% mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd)[, 3:4]) %>% mutate(transformation = "VarianceStabilizingTransformation"),
  as_tibble(assay(rld)[, 3:4]) %>% mutate(transformation = "rLogTransformation (DESeq2 default)"))

colnames(df)[1:2] <- c("x", "y")
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation) + labs(x="rep1.RED", y="rep2.RED")
invisible(dev.off())

# Run the DESeq pipeline
dds <- DESeq(dds)

# Plot dispersions
cairo_pdf("RNASEq.WHITE-RED.dispersions.pdf",width=15,height=10,antialias="default",fallback_resolution = 300)
plotDispEsts(dds, main="Dispersion plot")
invisible(dev.off())

# Colors for plots below
## Use RColorBrewer, better
suppressPackageStartupMessages(library(RColorBrewer))
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
suppressPackageStartupMessages(library(gplots))
cairo_pdf("RNASEq.WHITE-RED.heatmap.pdf",width=15,height=10,antialias="default",fallback_resolution = 300)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
invisible(dev.off())

# Principal components analysis with pcaExplorer
suppressPackageStartupMessages(require(pcaExplorer))

cairo_pdf("RNASeq.WHITE-RED.PCA.pdf",width=15,height=10,antialias="default",fallback_resolution = 300,onefile=T)
pcaplot(rld,intgroup="condition",title = "RlogStabilizingTransformation (default DESeq2 method)")
pcaplot(vsd,intgroup="condition",title = "VarianceStabilizingTransformation")
invisible(dev.off())

# Get differential expression results
res <- results(dds)

# assign ISO gene names to gene IDs
suppressPackageStartupMessages(library(plyranges))

# Load genome annotation. It should be the same as used in RNASeq alignments!
gr  <- read_gff("/usr/local/genomes/hg38.gtf")
gr  <- filter (gr, type == "gene" ) %>% select(gene_id, gene_name)
out <- as.data.frame(gr)

# only gene IDs and names
out <- out %>% select (gene_id, gene_name)
out <- out %>% mutate (gene_id = gsub("\\.\\d+$","",gene_id))

# replace NAs in gene_name to gene_ids
out <- out %>% mutate (gene_name = coalesce(gene_name, gene_id))
m <- match(out$gene_id, rownames(res))
res.sub <- res[m,]
res.sub$symbol <- out$gene_name
mcols(res.sub)[7,] <- DataFrame(type="GeneName",description="Common Gene Name")
res <- res.sub

table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]

## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "GeneID"

resdata <- resdata %>% select(GeneID, symbol, everything())
names(resdata)[2]<- "GeneName"

## Write results
write.table(resdata, file="WHITE-RED.results.tsv",sep='\t')

cairo_pdf("WHITE-RED.maplot.pdf",width=15,height=10,pointsize=16,antialias="default",fallback_resolution = 300,onefile=T)
DESeq2::plotMA(res, ylim=c(-1,1))
invisible(dev.off())

suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("EnhancedVolcano"))

cairo_pdf("RNASeq.WHITE-RED.volcanoplot.pdf",width=15,height=10,pointsize=16,antialias="default",fallback_resolution = 300,onefile=T)

EnhancedVolcano(res,lab=res$symbol,x='log2FoldChange',y='pvalue',
    xlab = bquote(~Log[2]~ 'fold change'),
    pCutoff = 0.05,
    FCcutoff = 1.0,
#    pointSize = 0.5,
    labSize = 4.0,
    colAlpha = 1,
    legendLabels=c('NS','Log (base 2) fold-change','P value', 'P value & Log (base 2) fold-change'),
    legendPosition = 'bottom',
    legendLabSize = 10,
    legendIconSize = 3.0
#    drawConnectors = TRUE,
#    widthConnectors = 0.2,
#    colConnectors = 'grey30',
) + coord_cartesian(xlim=c(-5.5,5.5),ylim=c(0,78)) + scale_x_continuous(breaks=seq(-15,30,5), minor_breaks=seq(-14,27,1)) 

EnhancedVolcano(res,lab=res$symbol,x='log2FoldChange',y='padj',
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~-Log[10]~adjusted~italic(P)),
    pCutoff = 0.05,
    FCcutoff = 1.0,
#    pointSize = 0.5,
    labSize = 4.0,
    colAlpha = 1,
    legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
    legendPosition = 'bottom',
    legendLabSize = 10,
    legendIconSize = 3.0
#    drawConnectors = TRUE,
#    widthConnectors = 0.2,
#    colConnectors = 'grey30'
) + coord_cartesian(xlim=c(-5.5,5.5),ylim=c(0,70)) + scale_x_continuous(breaks=seq(-15,30,5), minor_breaks=seq(-14,27,1)) 

invisible(dev.off())
