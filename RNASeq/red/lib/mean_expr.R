#!/usr/bin/Rscript
# script to obtain mean expression values in TPM and FPKM from RSEM calculated tables.
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Dependency tools:
# R with libraries     tibble, dplyr
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
namehead = "K562.red.RNAseq"
nametail = ".genes.results"
# replicates number
REPs     = 2

for (num in 1:REPs) {
   fname <- paste0(namehead,num,nametail)
   temptbl <- read.table(fname, skip=1, sep = "\t", strip.white = TRUE)
   newname <- paste0("Expr",num)
   FPKM    <- paste0("FPKM",num)
   colnames(temptbl)[6] <- newname
   colnames(temptbl)[7] <- FPKM
   if (num == 1) {
         alldata <- as.data.frame (temptbl$newname,row.names=temptbl$V1)
         alldata <- as.data.frame (temptbl[[newname]],row.names=as.vector(temptbl$V1))
         colnames(alldata)[num] <- newname

         FPdata <- as.data.frame (temptbl$FPKM,row.names=temptbl$V1)
         FPdata <- as.data.frame (temptbl[[FPKM]],row.names=as.vector(temptbl$V1))
         colnames(FPdata)[num] <- FPKM
         next
                 }
   alldata <- add_column(alldata, temptbl[[newname]])
   FPdata  <- add_column(FPdata,  temptbl[[FPKM]])
   colnames(alldata)[num] <- newname
   colnames(FPdata)[num]  <- FPKM
                 }

  medvect  <- apply(alldata, 1, median)
#  meanvect <- apply(alldata, 1, mean, trim=0.2) # trimmed mean
  meanvect <- apply(alldata, 1, mean)          # "true" mean

meanpkm  <- apply(FPdata, 1, mean)          # FPKM "true" mean
medpkm   <- apply(FPdata, 1, median)        # FPKM median

alldata <- add_column(alldata,medvect,meanvect)
alldata <- rownames_to_column(alldata, var = "GeneID") %>% as_tibble()
finalmed  <- alldata %>% select(GeneID,medvect)
finalmean <- alldata %>% select(GeneID,meanvect)

FPdata <-  add_column(FPdata,medpkm,meanpkm)
FPdata <-  rownames_to_column(FPdata, var = "GeneID") %>% as_tibble()
finalPKmed  <- FPdata %>% select(GeneID,medpkm)
finalPKmean <- FPdata %>% select(GeneID,meanpkm)

colnames(finalPKmean)[2] <- "FPKM, mean"
colnames(finalmean)[2]   <- "TPM, mean"

write.table(as.data.frame(alldata), file="K562.red.hg38.TPM", row.names=F, col.names=T, sep="\t", quote=F)
write.table(as.data.frame(FPdata), file="K562.red.hg38.FPKM", row.names=F, col.names=T, sep="\t", quote=F)
write.table(as.data.frame(finalmean), file="K562.red.hg38.mean.TPM", row.names=F, col.names=T, sep="\t", quote=F)
write.table(as.data.frame(finalPKmean), file="K562.red.hg38.mean.FPKM", row.names=F, col.names=T, sep="\t", quote=F)
