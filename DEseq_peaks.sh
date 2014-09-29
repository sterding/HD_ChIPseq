#######################################
# Script to study differential H3k4me3 signal at enhancer regions, using DESeq2
# Author: Xianjun Dong
# Date: 2013-12-14
# Version: 0.0
#######################################

#!/bin/sh

cd ~/projects/HD/data/scratch

echo "Step 1: calculate reads count of h3k4me3 in the peak regions"
#######################################################################

for i in h3k4me3_HD_batch1_*/accepted_hits.bam Schraham_*/accepted_hits.bam; do
    qsub ~/projects/HD/src/_get_readscount_per_gene.sge $i $enhancers $i.RCinpeaks
done

## wait until all jobs are done

echo "Step 2: cat the readscount file together"
#######################################################################

RC_in_enhancer=~/projects/HD/data/h3k4me3.in.enhancer.allsamples.txt
cut -f1-6 ~/scratch/combined_Lane1_4031_1/*RCinEnhancer > tmp
echo -en "chr\tstart\tend\tID\tlength\tstrand" > $RC_in_enhancer
for i in ~/scratch/combined_Lane*/*RCinEnhancer ~/scratch/h3k4me3_HD_batch1*/*RCinEnhancer ~/scratch/Schraham_*/*RCinEnhancer;
do
    echo -ne "\t"`echo $i | sed 's/.*scratch\///;s/\/.*//'` >> $RC_in_enhancer
    cut -f7 $i | paste tmp - > tmp0
    mv tmp0 tmp
done
echo >> $RC_in_enhancer
cat tmp >> $RC_in_enhancer
#rm tmp

echo "Step 3: run DESeq2 to calculate the DE genes"
#######################################################################

## R script

countData <- read.table('~/projects/HD/data/h3k4me3.in.enhancer.allsamples.txt', header=T)
rownames(countData)=countData[,4]
countData=countData[,-c(1:6)]
countData=countData[,grep("h3k4me3|R7|3706R|7244|R30|1644", names(countData))] # only take the 6 batch1 and 5 cases of Male and g43yr
colData <- data.frame(condition=c(rep("case", sum(grepl("combined|h3k4me3", names(countData)))), rep("control", sum(grepl("Schraham_", names(countData))))))
rownames(colData)=names(countData)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
colData(dds)$condition <- factor(colData(dds)$condition, levels=c("control","case"))
colData(dds)$condition <- relevel(colData(dds)$condition, "control")

dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]
#head(res)

write.table(as.data.frame(res), file="../result/h3k4me3.in.enhancer.allsamples.DESeq2.xls", quote =F, sep="\t", col.names = NA, row.names = TRUE)

# plot
pdf("../result/h3k4me3.in.enhancer.allsamples.DESeq2.pdf", width=6, height=6)
plotMA(dds,ylim=c(-2,2),main="DESeq2", pvalCutoff=.0001, pointcol = c("#00000020","#ff0000"))
legend('topleft', "adjusted p<0.0001", col='red', pch=20, bty='n')
#plot(res$log2FoldChange, -log10(res$padj), pch=20, cex=.5,  col='#99999966')
dev.off()

echo "Step 4: GO analysis of DE genes"
#######################################################################
awk '$7<0.0001' ../result/h3k4me3.in.enhancer.allsamples.DESeq2.xls | cut -f1 | sed 's/___/\t/g' | cut -f2-4 | sort -u > ../result/h3k4me3.in.enhancer.allsamples.DESeq2.p0.0001.txt