### [optional] to check the difference btw v1 and v2 expression data from Adam
df=read.table("allgenes.H3k4me3.expression.v1")
rownames(df)=df[,1]; df=df[,-c(1,8)]; #,8,17:19,)];
colnames(df)=c("chr", "start_promoter", "end_promoter", "transID", "type","name", "start_peak", "end_peak", "length", "summit_pos", "tags_in_peak","-10log10(pvalue)","fold_enrichment", "FDR(%)", "summit_height", "peak_type", "raw_read_1", "raw_read_2", "M","A", "pvalue_MAnorm", "dis_tss2sumit", "C_0002", "C_0003", "C_0004", "C_0005", "C_0006", "C_0008", "C_0009", "C_0010", "C_0011", "C_0012", "C_0013", "C_0014_HD", "C_0014_PD", "C_0015", "C_0016", "C_0017", "C_0018_HD", "C_0018_PD", "C_0020", "C_0021_HD", "C_0021_PD", "C_0022", "C_0023", "C_0024", "C_0025", "C_0026", "C_0029_HD", "C_0029_PD", "C_0031", "C_0032", "C_0033", "C_0035", "C_0036", "C_0037", "C_0038", "C_0039", "C_0050", "C_0053", "C_0060", "C_0061", "C_0062", "C_0065", "C_0069", "C_0070", "C_0071", "C_0075", "C_0076", "C_0077", "C_0081", "C_0082", "C_0083", "C_0087", "H_0001", "H_0002", "H_0003", "H_0005", "H_0006", "H_0007", "H_0008", "H_0009", "H_0010", "H_0012", "H_0013", "H_0539", "H_0657", "H_0658", "H_0681", "H_0695", "H_0700", "H_0726", "H_0740", "H_0750");
# remove _HD control
df=subset(df, chr!="chrM", select=grep("_HD", colnames(df), invert=T))
#write.table(df, "~/projects/bu_neuro/data/allgenes.H3k4me3.expression.v1.tab",  quote =F, sep="\t", col.names = NA, row.names = TRUE)
df1=df

df=read.table("allgenes.H3k4me3.expression.v2")
rownames(df)=df[,1]; df=df[,-c(1,8)]; #,8,17:19,)];
colnames(df)=c("chr", "start_promoter", "end_promoter", "transID", "type","name", "start_peak", "end_peak", "length", "summit_pos", "tags_in_peak","-10log10(pvalue)","fold_enrichment", "FDR(%)", "summit_height", "peak_type", "raw_read_1", "raw_read_2", "M","A", "pvalue_MAnorm", "dis_tss2sumit", "C_0002", "C_0003", "C_0004", "C_0005", "C_0006", "C_0008", "C_0009", "C_0010", "C_0011", "C_0012", "C_0013", "C_0014_PD", "C_0015", "C_0016", "C_0017", "C_0018_PD", "C_0020", "C_0021_PD", "C_0022", "C_0023", "C_0024", "C_0025", "C_0026", "C_0029_PD", "C_0031", "C_0032", "C_0033", "C_0035", "C_0036", "C_0037", "C_0038", "C_0039", "C_0050", "C_0053", "C_0060", "C_0061", "C_0062", "C_0065", "C_0069", "C_0070", "C_0071", "C_0075", "C_0076", "C_0077", "C_0081", "C_0082", "C_0083", "C_0087", "H_0001", "H_0002", "H_0003", "H_0005", "H_0006", "H_0007", "H_0008", "H_0009", "H_0010", "H_0012", "H_0013", "H_0539", "H_0657", "H_0658", "H_0681", "H_0695", "H_0700", "H_0726", "H_0740", "H_0750");
df=subset(df, chr!="chrM", select=grep("_HD", colnames(df), invert=T))
df2=df
pdf("v1_v2.pdf", , colormodel='cmyk')
for(i in colnames(df1)[grep("^[H|C]_", colnames(df1))][40:68]){
    #print(i)
    png(paste("../results/h3k4me3/expression_v1v2",i,"png",sep="."))
    plot(df1[,i]+1, df2[,i]+1, asp=1, xlab="before", ylab="after", pch=20, cex=.5, col="#000000", log="xy", main=i)
    dev.off()
}
dev.off()

## new version

df=read.table("allgenes.H3k4me3inpromoter.expression.tab")
rownames(df)=df[,1]; df=df[,-c(1,8:10,19)];
colnames(df)=c("chr", "start_promoter", "end_promoter", "transID", "type","name", "start_peak", "end_peak", "peak_type", "raw_read_1", "raw_read_2", "M","A", "pvalue_MAnorm", "dis_tss2peakcentral", "C_0002", "C_0003", "C_0004", "C_0005", "C_0006", "C_0008", "C_0009", "C_0010", "C_0011", "C_0012", "C_0013", "C_0014_PD", "C_0015", "C_0016", "C_0017", "C_0018_PD", "C_0020", "C_0021_PD", "C_0022", "C_0023", "C_0024", "C_0025", "C_0026", "C_0029_PD", "C_0031", "C_0032", "C_0033", "C_0035", "C_0036", "C_0037", "C_0038", "C_0039", "C_0050", "C_0053", "C_0060", "C_0061", "C_0062", "C_0065", "C_0069", "C_0070", "C_0071", "C_0075", "C_0076", "C_0077", "C_0081", "C_0082", "C_0083", "C_0087", "H_0001", "H_0002", "H_0003", "H_0005", "H_0006", "H_0007", "H_0008", "H_0009", "H_0010", "H_0012", "H_0013", "H_0539", "H_0657", "H_0658", "H_0681", "H_0695", "H_0700", "H_0726", "H_0740", "H_0750");
df=subset(df, chr!="chrM", select=grep("_HD", colnames(df), invert=T))

## 1. expression of genes with h3k4me3 peaks vs. w/o h3k4me3 peak
# defining cutoff for lowly / highly expressed genes
pdf("highlow_expression.pdf", colormodel='cmyk', width=12, height=5)
par(mfrow=c(1,2))
hist(log2(rowMeans(df[,grep("H_", colnames(df))])), breaks=20, main='Distribution of mean expression of HD samples')
hist(log2(rowMeans(df[,grep("C_", colnames(df))])), breaks=20, main='Distribution of mean expression of control samples')
#par(mfrow=c(1,2))
#hist(rowMeans(log2(1+df[,grep("H_", colnames(df))])), breaks=20, main='Distribution of mean expression of HD samples')
#hist(rowMeans(log2(1+df[,grep("C_", colnames(df))])), breaks=20, main='Distribution of mean expression of control samples')
dev.off()
# chisq test
chisq.test(table(log2(rowMeans(df[,grep("H_", colnames(df))]))>5, df$start_peak==-1))
# X-squared = 19041.43, df = 1, p-value < 2.2e-16
chisq.test(table(log2(rowMeans(df[,grep("C_", colnames(df))]))>5, df$start_peak==-1))
# X-squared = 19553.78, df = 1, p-value < 2.2e-16

#df1=data.frame(HD_mean=log2(1+rowMeans(df[,grep("H_", colnames(df))])), H3K4me3.peak.in.promoter=ifelse(df$start_peak==-1, 'without.H3K4me3.peak.in.promoter','with.H3K4me3.peak.in.promoter'))
#df1=data.frame(HD_mean=rowMeans(log2(1+df[,grep("H_", colnames(df))])), H3K4me3.peak.in.promoter=ifelse(df$start_peak==-1, 'w/o H3K4me3 peak\nin promoter','w/ H3K4me3 peak\nin promoter'))
pdf("highlow_expression_h3k4me3.pdf")
df1=data.frame(HD_mean=rowMeans(log2(1+df[,grep("H_", colnames(df))])), H3K4me3.peak.in.promoter=grepl("peak1|common", df$peak_type))
boxplot(HD_mean~H3K4me3.peak.in.promoter, df1, ylab="Mean expression of HD samples", xaxt="n")
axis(1, at=c(1,2),adj=1,padj=0.5,labels=c('w/o H3K4me3 peak\nin promoter','w/ H3K4me3 peak\nin promoter'))
dev.off()

## 2. differnetially expressed gene with peak vs. w/o peak
df22=log2(1+df[,grep("^[H|C]_0", colnames(df))])
df2=df22[rowMeans(df22)!=0,]  # 39940 out of 57244 remained (69.7%)
defactor <- function(x) {return(as.numeric(levels(x))[x]);}
#summit_height=defactor(df$A[rowMeans(df22)!=0])
A=defactor(df$A[rowMeans(df22)!=0])
M=defactor(df$M[rowMeans(df22)!=0])
r1=M/2+A #df$raw_read_1[rowMeans(df22)!=0]
r2=A-M/2 #df$raw_read_2[rowMeans(df22)!=0]
peak_type=as.character(df$peak_type[rowMeans(df22)!=0])
summit_height=ifelse(peak_type=="unique_peak1", r1, ifelse(peak_type=="unique_peak2", r2, ifelse(grepl("common",peak_type), A, NA)))

#FC=apply(df,1,function(x) {(mean(x[grep("H_", colnames(df))])+1) / (1+mean(x[grep("C_", colnames(df))]))})
#pvalues=-log10(apply(df,1,function(x) {t.test(x[grep("H_", colnames(df))], x[grep("C_", colnames(df))])$p.value}))

# sum(apply(df2, 1, function(x) {mean(x[grep("C_", colnames(df2))])==mean(x[grep("H_", colnames(df2))])})) # none
pvalues=-log10(apply(df2,1,function(x) {t.test(x[grep("H_", colnames(df2))], x[grep("C_", colnames(df2))])$p.value}))
M=apply(df2,1,function(x) {mean(x[grep("H_", colnames(df2))]) - mean(x[grep("C_", colnames(df2))])})
df2=cbind(df[rowMeans(df22)!=0,], expression_FC=M, expression_pvalue=pvalues)

df20=subset(df2, select=c(peak_type, expression_FC, expression_pvalue), peak_type!=".")
#df20$peak_type=as.character(levels(df20$peak_type))[df20$peak_type]
df20$peak_type=gsub("common_peak[0-9]", "common_peak", df20$peak_type)
boxplot(expression_FC~peak_type, df20, outline = F)


df3=cbind(round(df2,2), summit_height=summit_height, peak_type=peak_type, pvalues_DE=pvalues,M_DE=M)
df3=subset(df3, peak_type!=".")
write.table(df3, "genes.withH3k4me3inpromoter.expressed.tab",  quote =F, sep="\t", col.names = NA, row.names = TRUE)


# volcano plot with different type of peaks
pdf("volcano_plot_allgenes.pdf", colormodel='cmyk')
plot(M[peak_type=="."], pvalues[peak_type=="."], pch=20, cex=.5, col='#99999966', xlab="M (log(HD)-log(Control))", ylab="-log10(p-value)", xlim=range(M), ylim=range(pvalues))
points(M[grep("common", peak_type)], pvalues[grep("common", peak_type)], pch=20, cex=(summit_height[grep("common", peak_type)]/4), col='#0000ff66') # commmon peaks
points(M[peak_type=="unique_peak2"], pvalues[peak_type=="unique_peak2"], pch=20, cex=(summit_height[peak_type=="unique_peak2"]/3+0.001), col='#00ff0088') # unique to control
points(M[peak_type=="unique_peak1"], pvalues[peak_type=="unique_peak1"], pch=20, cex=(summit_height[peak_type=="unique_peak1"]/3+0.001), col='#ff000088') # unique to HD
abline(h=2, lty=2, lwd=.7)
text(2.5,2, "p-value=0.01", adj=c(0,-0.4))
n=c(sum(peak_type=="unique_peak1"), sum(peak_type=="unique_peak2"), sum(grepl("common", peak_type)), sum(peak_type=="."))
lg=legend("topleft", paste(c("Genes with HD-specific H3K4me3 peak in promoter", "Genes with Control-specific H3K4me3 peak in promoter", "Genes with H3K4me3 peak in both HD and Control", "Genes with no H3K4me3 peak in promoter"), " (n=", n,")", sep=""), col=c('#ff0000', '#00ff00', '#0000ff', '#000000'), pch=20, pt.cex=c(2,2,2,1), bty='n')
legend(lg$rect$left,lg$rect$top-lg$rect$h, title="H3K4me3 peak height", title.adj=0.5, c(50,100,200,400), pt.cex=(c(50,100,200,400)/200)^2+0.3, col="#999999", pch=20,  bty="n", horiz=T)

plot(rowMeans(df2[pvalues>2 & peak_type==".",grep("H_", colnames(df2))]), rowMeans(df2[pvalues>2 & peak_type==".",grep("C_", colnames(df2))]), pch=20, cex=0.5, xlim=range(rowMeans(df2[pvalues>2,])), ylim=range(rowMeans(df2[pvalues>2,])), xlab="Mean expression for HD samples (log)", ylab="Mean expression for Control samples (log)")
points(rowMeans(df2[pvalues>2 & grepl("common", peak_type),grep("H_", colnames(df2))]), rowMeans(df2[pvalues>2 & grepl("common", peak_type),grep("C_", colnames(df2))]), pch=20, cex=(summit_height[pvalues>2 & grepl("common", peak_type)]/200)^2+0.3, col='#0000ff66')
points(rowMeans(df2[pvalues>2 & peak_type=="unique_peak2",grep("H_", colnames(df2))]), rowMeans(df2[pvalues>2 & peak_type=="unique_peak2",grep("C_", colnames(df2))]), pch=20, cex=(summit_height[pvalues>2 & peak_type=="unique_peak2"]/200)^2+0.3, col='#00ff0088')
points(rowMeans(df2[pvalues>2 & peak_type=="unique_peak1",grep("H_", colnames(df2))]), rowMeans(df2[pvalues>2 & peak_type=="unique_peak1",grep("C_", colnames(df2))]), pch=20, cex=(summit_height[pvalues>2 & peak_type=="unique_peak1"]/200)^2+0.3, col='#ff000088')
n=c(sum(pvalues>2 & peak_type=="unique_peak1"), sum(pvalues>2 & peak_type=="unique_peak2"), sum(pvalues>2 & grepl("common", peak_type)), sum(pvalues>2 & peak_type=="."))
lg=legend("topleft", paste(c("DE genes with HD-specific H3K4me3 peak in promoter", "DE genes with Control-specific H3K4me3 peak in promoter", "DE genes with H3K4me3 peak in both HD and Control", "DE genes with no H3K4me3 peak in promoter"), " (n=", n,")", sep=""), col=c('#ff0000', '#00ff00', '#0000ff', '#000000'), pch=20, pt.cex=c(2,2,2,1), bty='n')
legend(lg$rect$left,lg$rect$top-lg$rect$h, title="H3K4me3 peak height", title.adj=0.5, c(50,100,200,400), pt.cex=(c(50,100,200,400)/200)^2+0.3, col="#999999", pch=20,  bty="n", horiz=T)

dev.off()
