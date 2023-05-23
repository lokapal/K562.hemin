#!/bin/bash
#hg19 2864785220
#hg38 2913022398
TYPE=white
bamCoverage --effectiveGenomeSize 2913022398 --normalizeUsing RPKM --ignoreForNormalization chr14 --exactScaling -p 30 -b $TYPE.rep1.bam -o rep1.bw
bamCoverage --effectiveGenomeSize 2913022398 --normalizeUsing RPKM --ignoreForNormalization chr14 --exactScaling -p 30 -b $TYPE.rep2.bam -o rep2.bw
multiBigwigSummary bins -p 12 -b rep1.bw rep2.bw -o res12hs.npz
plotCorrelation -in res12hs.npz -c spearman --removeOutliers --skipZeros -p scatterplot -o 4C.K562.$TYPE.spearman.pdf --labels 4C.$TYPE.rep1 4C.$TYPE.rep2 --log1p
plotCorrelation -in res12hs.npz -c pearson  --removeOutliers --skipZeros -p scatterplot -o 4C.K562.$TYPE.pearson.pdf  --labels 4C.$TYPE.rep1 4C.$TYPE.rep2 --log1p
plotPCA -in res12hs.npz -o K562.white.4C.PCA.pdf -T "PCA of K562 white 4C read normalized counts"
