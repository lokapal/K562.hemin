#!/usr/bin/Rscript
# script to perform differential 4C analysis 
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Input:  counts.txt pre-computed file by featureCounts, filtered to 100+ contacts
# Output: 1. K562.4C.WHITE-RED.results.tsv          K562 white-red 4C Differential analysis tabbed file
#         2. K562.4C.WHITE-RED.scatterplots.pdf     Scatterplots to control replicates consistency
#         3. K562.4C.WHITE-RED.volcanoplot.pdf      VolcanoPlot to display the most prominent 4C down/upcontacted genes
#         4. K562.4C.WHITE-RED.PCA.pdf              PCA to control replicates consistency
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
countdata <- read.table("counts.selected", header=TRUE, row.names=1)

# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,6:ncol(countdata)]

# Remove .bam or .sam from filenames
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))
colnames(countdata) <- gsub("\\.\\.\\.", "", colnames(countdata))

# Convert to matrix
countdata <- as.matrix(countdata)

# Assign condition (first two are expansion, second two are control)
(condition <- factor(c(rep("ctl", 2), rep("exp", 2))))

# Analysis with DESeq2 ----------------------------------------------------

suppressPackageStartupMessages(library(DESeq2))

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)

vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
rld <- rlogTransformation(dds,fitType="local")
hist(assay(rld))

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))

ddc <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(ddc, normalized=TRUE)[, 1:2]+1)) %>% mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "VarianceStabilizingTransformation"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rLogTransformation (DESeq2 default)"))

colnames(df)[1:2] <- c("x", "y")  

cairo_pdf("K562.4C.scatterplots.pdf",width=15,height=10,antialias="default",fallback_resolution = 300,onefile=T)
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation) + labs(x = "WHITE.rep1", y="WHITE.rep2")

rm (df)
df <- bind_rows(
  as_data_frame(log2(counts(ddc, normalized=TRUE)[, 3:4]+1)) %>% mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 3:4]) %>% mutate(transformation = "VarianceStabilizingTransformation"),
  as_data_frame(assay(rld)[, 3:4]) %>% mutate(transformation = "rLogTransformation (DESeq2 default)"))

colnames(df)[1:2] <- c("x", "y")
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation) + labs(x="RED.rep1", y="RED.rep2")
invisible(dev.off())

# Run the DESeq pipeline, fitType="local" for 4Cseq data
dds <- DESeq(dds,fitType="local")

# Plot dispersions
cairo_pdf("K562.4C.qc-dispersions.pdf",width=15,height=10,antialias="default",fallback_resolution = 300)
plotDispEsts(dds, main="Dispersion plot")
invisible(dev.off())

## Use RColorBrewer, better
suppressPackageStartupMessages(library(RColorBrewer))
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
suppressPackageStartupMessages(library(gplots))
cairo_pdf("K562.4C.qc-heatmap.pdf",width=15,height=10,antialias="default",fallback_resolution = 300)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
invisible(dev.off())

suppressPackageStartupMessages(require(pcaExplorer))
cairo_pdf("K562.4C.PCA.pdf",width=15,height=10,antialias="default",fallback_resolution = 300,onefile=T)
pcaplot(rld,intgroup="condition",title = "RlogStabilizingTransformation 4C (default DESeq2 method)")
pcaplot(vsd,intgroup="condition",title = "VarianceStabilizingTransformation 4C")
invisible(dev.off())

# Get differential expression results
res <- results(dds)

# assign ISO gene names to gene IDs
suppressPackageStartupMessages(library(plyranges))

# Load genome annotation. It should be the same as used in 4C gene assignments!
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
write.table(resdata, file="K562.4C.WHITE-RED.results.tsv",sep='\t')

cairo_pdf("K562.4C.maplot.pdf",width=15,height=10,pointsize=16,antialias="default",fallback_resolution = 300,onefile=T)
plotMA(dds, ylim=c(-1,1), cex=1)
invisible(dev.off())

suppressPackageStartupMessages(library("EnhancedVolcano"))

cairo_pdf("K562.4C.volcanoplot.pdf",width=15,height=10,pointsize=16,antialias="default",fallback_resolution = 300,onefile=T)
EnhancedVolcano(res,lab=res$symbol,x='log2FoldChange',y='pvalue',
    xlab = bquote(~Log[2]~ 'fold change'),
    pCutoff = 0.01,
    FCcutoff = 2,
    transcriptLabSize = 4.0,
    colAlpha = 1,
    legend=c('NS','Log (base 2) fold-change','P value', 'P value & Log (base 2) fold-change'),
    legendPosition = 'bottom',
    legendLabSize = 10,
    legendIconSize = 3.0,
#    drawConnectors = TRUE,
#    widthConnectors = 0.2,
#    colConnectors = 'grey30',
) + coord_cartesian(xlim=c(-21,21),ylim=c(0,300)) + scale_x_continuous(breaks=seq(-25,30,5), minor_breaks=seq(-26,26,1)) + scale_y_continuous(breaks=seq(0,300,50),minor_breaks=seq(1,500,10))

EnhancedVolcano(res,lab=res$symbol,x='log2FoldChange',y='padj',
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~-Log[10]~adjusted~italic(P)),
    pCutoff = 0.01,
    FCcutoff = 2,
    transcriptLabSize = 4.0,
    colAlpha = 1,
    legend=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
    legendPosition = 'bottom',
    legendLabSize = 10,
    legendIconSize = 3.0,
#    drawConnectors = TRUE,
#    widthConnectors = 0.2,
#    colConnectors = 'grey30'
) + coord_cartesian(xlim=c(-21,21),ylim=c(0,300)) + scale_x_continuous(breaks=seq(-25,30,5), minor_breaks=seq(-26,26,1)) + scale_y_continuous(breaks=seq(0,300,50),minor_breaks=seq(1,500,10))

invisible(dev.off())
