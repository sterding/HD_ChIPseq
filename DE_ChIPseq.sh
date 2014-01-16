#######################################
# Script to study differential H3k4me3 signal at enhancer regions, using DESeq2
# Author: Xianjun Dong
# Date: 2013-12-14
# Version: 0.0
#######################################

#!/bin/sh

# required files
annotation=$GENOME/hg19/Annotation/Genes/gencode.v17.annotation.gtf.gz
enhancers=~/projects/bu_neuro/data/enhancer.superenhancer.hg19.BI_Brain_Mid_Frontal_Lobe.bed
#grep -v track ~/projects/PD/superEnhancer/enhancers.superenhancer.hg19.DataS2/BI_Brain_Mid_Frontal_Lobe.bed | sort -k5,5n -u | awk '{OFS="\t";print $0,"."}' > ~/projects/HD/data/enhancer.superenhancer.hg19.BI_Brain_Mid_Frontal_Lobe.bed

echo "Step 1: calculate reads count of h3k4me3 in the enhancer regions"
#######################################################################

for i in ~/scratch/combined_Lane*/accepted_hits.bam ~/scratch/h3k4me3_HD_batch1*/accepted_hits.bam ~/scratch/Schraham_*/accepted_hits.bam; do
    qsub ~/projects/bu_neuro/src/_get_readscount_per_gene.sge $i $enhancers $i.RCinEnhancer
done

## wait until all jobs are done

echo "Step 2: cat the readscount file together"
#######################################################################

RC_in_enhancer=~/projects/bu_neuro/data/h3k4me3.in.enhancer.allsamples.txt
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