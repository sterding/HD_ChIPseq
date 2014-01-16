#######################################
# Script to call gene/tx with differential H3k4me3 signal at promoter regions, using DESeq2
# Author: Xianjun Dong
# Date: 2013-12-14
# Version: 0.0
#######################################

#!/bin/sh

# required files
annotation=$GENOME/hg19/Annotation/Genes/gencode.v17.annotation.gtf.gz
promoters=$annotation.promoter_1k_round_TSS.bed
#zcat $annotation | fgrep -w transcript | sed 's/[";]//g;' | awk '{OFS="\t"; tss=($7=="+")?($4-1):($5-1); if(tss<2000) tss=2000; print $1, tss-2000,tss+2000,$12"___"$10"___"$14"___"$18, 4000, $7}'  > $annotation.promoter_4k_round_TSS.bed
#zcat $annotation | fgrep -w transcript | sed 's/[";]//g;' | awk '{OFS="\t"; tss=($7=="+")?($4-1):($5-1); if(tss<500) tss=500; print $1, tss-500,tss+500,$12"___"$10"___"$14"___"$18, 1000, $7}'  > $annotation.promoter_1k_round_TSS.bed

echo "Step 1: calculate reads count of h3k4me3 in the promoter regions"
#######################################################################

for i in ~/scratch/combined_Lane*/accepted_hits.bam ~/scratch/h3k4me3_HD_batch1*/accepted_hits.bam ~/scratch/Schraham_*/accepted_hits.bam; do
    qsub ~/projects/bu_neuro/src/_get_readscount_per_gene.sge $i $promoters $i.RCinPromoter
done

## wait until all jobs are done

echo "Step 2: cat the readscount file together"
#######################################################################

RC_in_promoter=~/projects/bu_neuro/data/h3k4me3.in.promoter_1k_round_TSS.allsamples.txt
cut -f1-6 ~/scratch/combined_Lane1_4031_1/*RCinPromoter > tmp
echo -en "chr\tstart\tend\tID\tlength\tstrand" > $RC_in_promoter
for i in ~/scratch/combined_Lane*/*RCinPromoter ~/scratch/h3k4me3_HD_batch1*/*RCinPromoter ~/scratch/Schraham_*/*RCinPromoter;
do
    echo -ne "\t"`echo $i | sed 's/.*scratch\///;s/\/.*//'` >> $RC_in_promoter
    cut -f7 $i | paste tmp - > tmp0
    mv tmp0 tmp
done
echo >> $RC_in_promoter
cat tmp >> $RC_in_promoter
#rm tmp

echo "Step 3: run DESeq2 to calculate the DE genes"
#######################################################################

## R script

countData <- read.table('~/projects/HD/data/h3k4me3.in.promoter_1k_round_TSS.allsamples.txt', header=T)
rownames(countData)=countData[,4]
countData=countData[,-c(1:6)]

# similarity btw samples
dd <- as.dist((1 - cor(countData))/2)
round(1000 * dd) # (prints more nicely)
plot(hclust(dd)) # to see a dendrogram of clustered variables

dd <- as.dist((1 - cor(countData[,grep("Schraham", colnames(countData), invert=T)]))/2)
round(1000 * dd) # (prints more nicely)
plot(hclust(dd)) # to see a dendrogram of clustered variables


countData=countData[,grep("h3k4me3|R7|3706R|7244|R30|1644", names(countData))] # only take the 6 batch1 and 5 cases of Male and g43yr
#colData <- data.frame(condition=c(rep("case", sum(grepl("combined|h3k4me3", names(countData)))), rep("control", sum(grepl("Schraham_", names(countData))))))
#rownames(colData)=names(countData)
# incl covairance
colData=read.table("../data/covariances.tab", header=T)
rownames(colData)=colData[,1]
colData=colData[,-1]

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition) # + age + gender + PMI)
colData(dds)$condition <- factor(colData(dds)$condition, levels=c("control","case"))
colData(dds)$condition <- relevel(colData(dds)$condition, "control")

dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]
#head(res)

write.table(as.data.frame(res), file="../result/h3k4me3.in.promoter_1k_round_TSS.allsamples.DESeq2.xls", quote =F, sep="\t", col.names = NA, row.names = TRUE)

# plot
png("../results/h3k4me3.in.promoter_1k_round_TSS.allsamples.DESeq2.png", width=800, height=800)
plotMA(dds,ylim=c(-2,2),main="DESeq2", pvalCutoff=.0001, pointcol = c("#00000020","#ff0000"))
legend('topleft', "genes with adjusted p<0.0001", col='red', pch=20, bty='n')
#res[is.na(res$pvalue),c('pvalue','padj')]=1
#res=subset(res, baseMean>0)
#plot(res$baseMean, pmax(-2, pmin(2, res$log2FoldChange)), log = 'x', pch = ifelse(res$log2FoldChange < -2, 6, ifelse(res$log2FoldChange > 2, 2, 20)), cex=.5, col=ifelse(res$padj<.0001, '#ff0000', '#00000020'))
#abline(h = 0, lwd = 4, col = '#ff000080')

## color with gene expression fold change!!!
# take the tx with the highest promoter signal as representative (done in bash)
# grep -p "^ENST" result/h3k4me3.in.promoter_1k_round_TSS.allsamples.DESeq2.xls | sed 's/___/ /g' | sort -k2,2 -k5,5nr | awk '!arr[$2]++' > result/h3k4me3.in.promoter_1k_round_TSS.allsamples.DESeq2.TX.xls
# merge with gene expression file (from RNAseq)
# 
res2 <- read.table('~/projects/HD/result/Dec16/h3k4me3.in.promoter_1k_round_TSS.allsamples.DESeq2.TX.xls', header=F)
names(res2)=c('Tx_ID', 'Gene_ID', 'Gene_Type', 'Gene_Symbol', "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
res2[is.na(res2$pvalue),'pvalue']=1
res2[is.na(res2$padj),'padj']=1
res2=subset(res2, baseMean>0)
rownames(res2)=res2$Gene_ID; res2=res2[,-2]

# read expression data
df=read.table("../data/allgenes.H3k4me3.expression.v2")
rownames(df)=df[,1]; df=df[,-c(1,8)]; #,8,17:19,)];
colnames(df)=c("chr", "start_promoter", "end_promoter", "transID", "type","name", "start_peak", "end_peak", "length", "summit_pos", "tags_in_peak","-10log10(pvalue)","fold_enrichment", "FDR(%)", "summit_height", "peak_type", "raw_read_1", "raw_read_2", "M","A", "pvalue_MAnorm", "dis_tss2sumit", "C_0002", "C_0003", "C_0004", "C_0005", "C_0006", "C_0008", "C_0009", "C_0010", "C_0011", "C_0012", "C_0013", "C_0014_PD", "C_0015", "C_0016", "C_0017", "C_0018_PD", "C_0020", "C_0021_PD", "C_0022", "C_0023", "C_0024", "C_0025", "C_0026", "C_0029_PD", "C_0031", "C_0032", "C_0033", "C_0035", "C_0036", "C_0037", "C_0038", "C_0039", "C_0050", "C_0053", "C_0060", "C_0061", "C_0062", "C_0065", "C_0069", "C_0070", "C_0071", "C_0075", "C_0076", "C_0077", "C_0081", "C_0082", "C_0083", "C_0087", "H_0001", "H_0002", "H_0003", "H_0005", "H_0006", "H_0007", "H_0008", "H_0009", "H_0010", "H_0012", "H_0013", "H_0539", "H_0657", "H_0658", "H_0681", "H_0695", "H_0700", "H_0726", "H_0740", "H_0750");
df=subset(df, chr!="chrM", select=grep("_HD", colnames(df), invert=T))
## TODO: ask Adam to provide statistics for the DEseq analysis, incl. p-value etc.

cn=intersect(rownames(df), rownames(res2))
res2=cbind(res2[cn,], mean_C=rowMeans(df[cn,grep("C_", colnames(df))]), mean_H=rowMeans(df[cn,grep("H_", colnames(df))]))
res2=cbind(res2, log2FC_expression=log2(.1+res2$mean_H)-log2(res2$mean_C+.1))
res2=cbind(res2, C=ifelse(log2(res2$mean_C+.1)>5, 1, -1), H=ifelse(log2(.1+res2$mean_H)>5, 1, -1))
res2=cbind(res2, HCmean=log2((res2$mean_C+res2$mean_H)/2))

h1=hist(log2(res2$baseMean[res2$HCmean<5]), breaks=seq(-5,10,.5));
h2=hist(log2(res2$baseMean[res2$HCmean>=5]), breaks=seq(-5,10,.5));
d=rbind(h1$count, h2$count)
colnames(d)=h1$mid
rownames(d)=c('low-expression', 'high-expression')

pdf("../result/2014Jan16/expression.vs.h3k4me3.pdf", width=7, height=4)
hist(log2(res2$mean_H), breaks=50, main='Distribution of RNAseq expression of HD samples')
hist(log2(res2$mean_C), breaks=50, main='Distribution of RNAseq expression of control samples')
hist(res2$HCmean, breaks=50, main='Distribution of RNAseq expression of all samples')
hist(log2(res2$baseMean), breaks=seq(-5,10,.5), ylab="Count", xlab="H3K4me3 Signal (log2 based)", main="Distribution of H3K4me3 signal")
barplot(d, main="Distribution by H3K4me3 and RNAseq signal", cex.name=.7, xlab="H3K4me3 Signal (log2 based)", ylab="Count", col=c("darkblue","red"),legend = rownames(d))
dev.off()

# DE gene vs. DP h3k4me3
pdf("../result/2014Jan16/DE.vs.DP.pdf", width=7, height=5)

plot(res2$baseMean, pmax(-2, pmin(2, res2$log2FoldChange)), log = 'x', xlab="mean of normalized counts", ylab="log2 fold change", pch = ifelse(res2$log2FoldChange < -2, 25, ifelse(res2$log2FoldChange > 2, 24, ifelse(res2$padj<.0001, 21,20))), cex=.5, lwd=.5, col=ifelse(res2$padj>=.0001, "#00000020", '#ff0000'), bg='#00000020')
abline(h = 0, lwd = 4, col = '#ff000080')
legend('topleft', "genes with differential H3K4me3 signal (FDR<0.0001)", col='red', bg="#00000020", pch=21, bty='n', cex=.7)

plot(res2$baseMean, pmax(-2, pmin(2, res2$log2FoldChange)), log = 'x', xlab="mean of normalized counts", ylab="log2 fold change", pch = ifelse(res2$log2FoldChange < -2, 25, ifelse(res2$log2FoldChange > 2, 24, ifelse(res2$padj<.0001, 21,20))), cex=.5, lwd=.5, col=ifelse(res2$padj>=.0001, "#00000020", '#ff0000'), bg=ifelse(res2$padj>=.0001, "gray", ifelse(res2$mean_H>res2$mean_C, 'green', 'blue')))
abline(h = 0, lwd = 4, col = '#ff000080')
legend('topleft', c("genes with differential H3K4me3 signal (FDR<0.0001)", "RNAseq: HD >= Ct", "RNAseq: HD < Ct"), col=rep('red', 3), pt.bg=c("#00000020", 'green', 'blue'), pch=rep(21, 3), bty='n', cex=.7)

plot(res2$baseMean, pmax(-2, pmin(2, res2$log2FoldChange)), log = 'x', xlab="mean of normalized counts", ylab="log2 fold change", pch = ifelse(res2$log2FoldChange < -2, 25, ifelse(res2$log2FoldChange > 2, 24, 21)), cex=.5, lwd=.5, col=ifelse(res2$padj>=.001, "#00000020", '#ff0000'), bg=ifelse(res2$mean_H>res2$mean_C, '#00ff0080', '#0000ff80'))
abline(h = 0, lwd = 4, col = '#ff000080')
legend('topleft', c("genes with differential H3K4me3 signal (FDR<0.0001)", "RNAseq: HD >= Ct", "RNAseq: HD < Ct"), col=c('red', "#00000020", "#00000020"), pt.bg=c("#00000020", '#00ff0080', '#0000ff80'), pch=rep(21, 3), bty='n', cex=.7)

dev.off()


plot(res2$baseMean, pmax(-2, pmin(2, res2$log2FoldChange)), log = 'x', pch = ifelse(res2$log2FoldChange < -2, 6, ifelse(res2$log2FoldChange > 2, 2, 20)), cex=.5, col=ifelse(res2$padj>=.0001, "#00000020", ifelse(res2$mean_H>res2$mean_C, 'green', 'blue')))
abline(h = 0, lwd = 4, col = '#ff000080')
legend('topleft', "genes with adjusted p<0.0001", col='red', pch=20, bty='n')

plot(res2$baseMean, pmax(-2, pmin(2, res2$log2FoldChange)), log = 'x', pch = ifelse(res2$log2FoldChange < -2, 6, ifelse(res2$log2FoldChange > 2, 2, ifelse(res2$padj<.0001, 21,20))), cex=.5, col=ifelse(res2$padj<.0001, '#ff0000', '#00000020'), bg=ifelse(abs(res2$log2FC_expression)<1, 'gray', ifelse(res2$log2FC_expression>1, 'green', 'blue')))

plot(res2$baseMean, pmax(-2, pmin(2, res2$log2FoldChange)), log = 'x', pch = ifelse(res2$log2FoldChange < -2, 6, ifelse(res2$log2FoldChange > 2, 2, ifelse(res2$padj<.0001, 21,20))), cex=.5, col=ifelse(res2$padj<.0001, '#ff0000', '#00000020'), bg=ifelse(abs(res2$C+res2$H)>0, 'gray',ifelse(res2$H>res2$C, 'green', 'blue')))

abline(h = 0, lwd = 4, col = '#ff000080')

dec.off()



echo "Step 4: GO analysis of DE genes"
#######################################################################
awk '$7<0.0001' ../result/h3k4me3.in.promoter_1k_round_TSS.allsamples.DESeq2.xls | cut -f1 | sed 's/___/\t/g' | cut -f2-4 | sort -u > ../result/h3k4me3.in.promoter_1k_round_TSS.allsamples.DESeq2.p0.0001.txt