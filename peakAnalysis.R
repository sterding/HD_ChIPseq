peaks=read.table("~/projects/HD/data/scratch/h3k4me3_peaks_union.signal.neighborhood.annotated.bed.xls", header=T)
attach(peaks)

x=log2(rpm_summit_CT+0.05)
y=log2(rpm_summit_HD+0.05)
# normalize
library('affy')
tmp = normalize.invariantset(x, y)
x.norm <- as.numeric(approx(tmp$n.curve$y, tmp$n.curve$x, xout = x, rule = 2)$y)  #https://support.bioconductor.org/p/3063/

# Question: how to define uniqueness?
uniqueness=ifelse(abs(x.norm-y)<=1, "HDCT", ifelse(y>x.norm,"HD","CT"))
proximality=ifelse(abs(distance_peakcenter2TSS)<=2000, "proximal", "distal")

#> table(uniqueness,proximality)
#          proximality
#uniqueness distal proximal
#      CT      424      462
#      HD      528      722
#      HDCT   7475    18997
     
# output to file for vann diagram
peaks=cbind(id=paste(chr,start_peak,end_peak, sep="_"), peaks, HD.lognorm=y, CT.lognorm=x.norm, uniqueness=uniqueness, proximality=proximality)
write.table(peaks, "~/Dropbox/huntington_bu/data/h3k4me3_peaks_union.signal.neighborhood.annotated.bed.v3.xls", quote =F, sep="\t", row.names =F, col.names =T)

# size distribution of peaks
pdf("~/Dropbox/huntington_bu/figs/peak_dize_distribution.pdf", width=6, height=4)
hist(pmin(5000,end_peak-start_peak), breaks=100)
hist(log10(pmin(5000,end_peak-start_peak)), breaks=100)
dev.off()


# distance of peak to TSS
####################################################
pdf("~/Dropbox/huntington_bu/figs/distance_peakcenter2TSS.pdf", width=6, height=4)
hist(pmax(pmin(distance_peakcenter2TSS, 1e4), -1e4), breaks=100, xlab="distance between the nearest TSS to the center of peak (bp)", main="histogram of distance from peak center to the nearest TSS")
hist(distance_peakcenter2TSS[distance_peakcenter2TSS < 2e3 & distance_peakcenter2TSS > -2e3], breaks=100, main="insert",xlab="")
abline(v=c(-1000,1000), col='red', lty=2)
dev.off()

#sample-specific peaks vs. universal peaks
pdf("~/Dropbox/huntington_bu/figs/peak_vs_sampleCount.pdf", width=6, height=6)
require('ggplot2')
df=data.frame(rpm_summit=0.1+c(rpm_summit_CT, rpm_summit_HD), samples_support_peak=c(samples_support_peak_CT, samples_support_peak_HD), group=rep(c("Healthy control", "Huntington's disease"), each=length(samples_support_peak_HD)))
p <- ggplot(df, aes(factor(samples_support_peak), rpm_summit))
p + geom_boxplot(aes(fill = factor(group)), outlier.size=1, outlier.colour='darkgray') +
    scale_y_log10(breaks=c(.1,1,10,100), labels=c(.1,1,10,100)) +
    labs(x="Number of samples supporting a peak", y='Normalized signal of the peak summit') +
    scale_fill_manual("Condition",values=c("#99d8c9","#fc9272")) +
    theme(legend.justification=c(0,1),legend.position=c(0,1), legend.background = element_rect(fill = 'transparent'), legend.key = element_rect(fill = 'transparent',colour = 'NA'))
plot(table(samples_support_peak_CT, samples_support_peak_HD), col=rev(gray.colors(7)), main="", xlab="number of CT samples supporting a peak", ylab="number of HD samples supporting a peak")
plot(table(samples_support_peak_CT, samples_support_peak_HD), col=1:7, off=0, main="", border='white', xlab="number of CT samples supporting a peak", ylab="number of HD samples supporting a peak")

dev.off()

# the following figure is to show that using cutoff of 3 samples to define unique peaks is as good as using a foldchange cutoff, except some outliers.
pdf("~/Dropbox/huntington_bu/figs/uniqness.foldchange.samplecount.pdf", width=7, height=6)
HD=rpm_summit_HD; CT=rpm_summit_CT; 
fold=ifelse(HD/CT>1,HD/CT,CT/HD)
ave=log((HD + CT)/2 + 1)
#ColorRamp=colorRampPalette(c("blue", "red"))(10)
col=rgb(ifelse(HD/CT>1,255,0),0,ifelse(HD/CT>1,0,255),250*(ave-min(ave))/(max(ave)-min(ave)), maxColorValue=255)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(samples_support_peak_CT, samples_support_peak_HD, cex=fold, col=col, main="no pseudocount", xlab="Number of CT samples called a peak", ylab="Number of HD samples called a peak", xlim=c(-1,6.7), ylim=c(-1,6.7),xaxp=c(0,6,6),yaxp=c(0,6,6))
legend("topright", inset=c(-0.3,0), legend=c(1:3)*5, pt.cex=c(1:3)*5, pch=1, bty='n', title="fold change")
x=c(0.5,1,5,10,100)
px=pmin(250,250*(log(x+1)-min(ave))/(max(ave)-min(ave)))
legend("bottomright", inset=c(-0.25,0), legend=c(paste(x,"rpm"), "HD>=CT","CT>HD"), pch=19, bty='n', col=c(rgb(0,0,0,px,maxColorValue=255),"red","blue"), title="signal value")

# a bigger pseudocount
HD=rpm_summit_HD+0.05; CT=rpm_summit_CT+0.05;  # 0.05 is a good number because it won't screw up the bio-model of log(summit_HD)
fold=ifelse(HD/CT>1,HD/CT,CT/HD)
ave=log((HD + CT)/2 + 1)
#ColorRamp=colorRampPalette(c("blue", "red"))(10)
col=rgb(ifelse(HD/CT>1,255,0),0,ifelse(HD/CT>1,0,255),250*(ave-min(ave))/(max(ave)-min(ave)), maxColorValue=255)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(samples_support_peak_CT, samples_support_peak_HD, cex=fold, col=col, main="pseudocount=0.05",xlab="Number of CT samples called a peak", ylab="Number of HD samples called a peak", xlim=c(-1,6.7), ylim=c(-1,6.7),xaxp=c(0,6,6),yaxp=c(0,6,6))
legend("topright", inset=c(-0.3,0), legend=c(1:3)*5, pt.cex=c(1:3)*5, pch=1, bty='n', title="fold change")
x=c(0.5,1,5,10,100)
px=pmin(250,250*(log(x+1)-min(ave))/(max(ave)-min(ave)))
legend("bottomright", inset=c(-0.25,0), legend=c(paste(x,"rpm"), "HD>=CT","CT>HD"), pch=19, bty='n', col=c(rgb(0,0,0,px,maxColorValue=255),"red","blue"), title="signal value")

dev.off()

table(fold>2, samples_support_peak_CT>=3 & samples_support_peak_HD>=3)

## uniq/common vs. proximal/distal
#1. 2x2 table


# differential peak btw HD and control
####################################################

# define unique/common: (a) >=2 fold between normalized signal or (b) >=3 difference in numbers of samples supporting the peak (e.g. 0|3, 1|4, 2|5, 3|6)  OR a simple cutoff (e.g. <3 vs. >=3)
table(foldChange=abs(x.norm-y)>=1,  sampleDiff=abs(samples_support_peak_CT-samples_support_peak_HD)>=3)
#          sampleDiff
#foldChange FALSE  TRUE
#     FALSE 24161  2311
#     TRUE   1214   922
# chisq.test p-value < 2.2e-16


pdf("~/Dropbox/huntington_bu/figs/peak_signal_HD.vs.Ct.pdf", width=7, height=7)

## XY plot
par(mfrow=c(2,2))
plot(x,y, pch='.', cex=1.5, col=ifelse(abs(distance_peakcenter2TSS)<=1000, "#3333ff22", "#ff333322"), ylim=range(x, y), xlim=range(x, y), xlab="Peak summit in control - log2(RPM)", ylab="Peak summit in HD - log2(RPM)", main="before normalization")
abline(a=0, b=1,  lty=2, lwd=1)
legend("topleft", c("proximal peaks", "distal peaks"), pch=20, col=c("#3333ff", "#ff3333"),  bty = "n")
plot(x.norm,y, pch='.', cex=1.5, col=ifelse(abs(distance_peakcenter2TSS)<=1000, "#3333ff22", "#ff333322"), ylim=range(x.norm, y), xlim=range(x.norm, y), xlab="Peak summit in control - log2(RPM)", ylab="Peak summit in HD - log2(RPM)", main="after normalization")
abline(a=0, b=1,  lty=2, lwd=1)
legend("topleft", c("proximal peaks", "distal peaks"), pch=20, col=c("#3333ff", "#ff3333"),  bty = "n")

## MAplot
M=y-x; A=(x+y)/2
plot(A,M, col=ifelse(abs(distance_peakcenter2TSS)<=1000, "#3333ff22", "#ff333322"), ylim=range(-max(M), max(M)), xlab="Mean log2 density (HD+Ct)", ylab="log2 foldchange (HD/Ct)", pch=".", cex=2, main="before normalization")
abline(h=0, col='red',lty=2, lwd=1)
M=y-x.norm; A=(x.norm+y)/2
plot(A,M, col=ifelse(abs(distance_peakcenter2TSS)<=1000, "#3333ff22", "#ff333322"), ylim=range(-max(M), max(M)), xlab="Mean log2 density (HD+Ct)", ylab="log2 foldchange (HD/Ct)", pch=".", cex=2, main="after normalization")
points(A[abs(M)>1 & A>0],M[abs(M)>1 & A>0], pch='.',cex=2, col=ifelse(abs(distance_peakcenter2TSS[abs(M)>1 & A>0])<=1000, "#0000ff", "#ff0000"))
abline(h=0, col='red',lty=2, lwd=1)

# separate into proximal and distal
# proximal
x=log2(rpm_summit_CT[abs(distance_peakcenter2TSS)<=1000]+0.05)
y=log2(rpm_summit_HD[abs(distance_peakcenter2TSS)<=1000]+0.05)
tmp = normalize.invariantset(x, y)
x.norm <- as.numeric(approx(tmp$n.curve$y, tmp$n.curve$x, xout = x, rule = 2)$y)  #https://support.bioconductor.org/p/3063/

# XYplot
plot(x,y, pch='.', cex=1.5, ylim=range(x, y), xlim=range(x, y), xlab="Peak summit in control - log2(RPM)", ylab="Peak summit in HD - log2(RPM)", main="before normalization")
abline(a=0, b=1, col='red', lty=2, lwd=1)
plot(x.norm,y, pch='.', cex=1.5, ylim=range(x.norm, y), xlim=range(x.norm, y), xlab="Peak summit in control - log2(RPM)", ylab="Peak summit in HD - log2(RPM)", main="after normalization")
abline(a=0, b=1, col='red', lty=2, lwd=1)

## MAplot
M=y-x; A=(x+y)/2
plot(A,M, ylim=range(-max(M), max(M)), xlab="Mean log2 density (HD+Ct)", ylab="log2 foldchange (HD/Ct)", pch=".", cex=1.5, main="before normalization")
abline(h=0, col='red',lty=2, lwd=1)
M=y-x.norm; A=(x.norm+y)/2
plot(A,M, ylim=range(-max(M), max(M)), xlab="Mean log2 density (HD+Ct)", ylab="log2 foldchange (HD/Ct)", pch=".", cex=1.5, main="after normalization")
points(A[abs(M)>1 & A>0],M[abs(M)>1 & A>0], pch=1, col='darkgreen')
abline(h=0, col='red',lty=2, lwd=1)

# distal
x=log2(rpm_summit_CT[abs(distance_peakcenter2TSS)>1000]+0.05)
y=log2(rpm_summit_HD[abs(distance_peakcenter2TSS)>1000]+0.05)
tmp = normalize.invariantset(x, y)
x.norm <- as.numeric(approx(tmp$n.curve$y, tmp$n.curve$x, xout = x, rule = 2)$y)  #https://support.bioconductor.org/p/3063/

plot(x,y, pch='.',cex=1.5,  ylim=range(x, y), xlim=range(x, y), xlab="Peak summit in control - log2(RPM)", ylab="Peak summit in HD - log2(RPM)", main="before normalization")
abline(a=0, b=1, col='red', lty=2, lwd=1)
plot(x.norm,y, pch='.',cex=1.5,  ylim=range(x.norm, y), xlim=range(x.norm, y), xlab="Peak summit in control - log2(RPM)", ylab="Peak summit in HD - log2(RPM)", main="after normalization")
abline(a=0, b=1, col='red', lty=2, lwd=1)

## MAplot
M=y-x; A=(x+y)/2
plot(A,M, ylim=range(-max(M), max(M)), xlab="Mean log2 density (HD+Ct)", ylab="log2 foldchange (HD/Ct)", pch=".", cex=1.5, main="before normalization")
abline(h=0, col='red',lty=2, lwd=1)
M=y-x.norm; A=(x.norm+y)/2
plot(A,M, ylim=range(-max(M), max(M)), xlab="Mean log2 density (HD+Ct)", ylab="log2 foldchange (HD/Ct)", pch=".", cex=1.5, main="after normalization")
points(A[abs(M)>1 & A>0],M[abs(M)>1 & A>0], pch=1, col='darkgreen')
abline(h=0, col='red',lty=2, lwd=1)

dev.off()

#
####################################################

q('no')

#################################################### (OLD VERSION) #############################################
peaks=read.table("significantPeaks.bed.inpromoter.tab", header=T)
#colnames(peaks)=c('chr',"start_peak", "end_peak", "length", "summit_pos", "tags_in_peak","pvalue_MACS","fold_enrichment", "FDR","summit_height","peak_type", "raw_read_1", "raw_read_2", "M","A", "pvalue_MAnorm","inpromoter")

pdf("h3k4me3_peaks.pdf", width=8, height=6, colormodel='cmyk')
par(mar=c(4,4,4,1))
layout(matrix(seq(2),nrow=1,ncol=2),widths=c(3,0.5),heights=1)
plot(peaks$A_rescaled, peaks$M_rescaled, col=rev(rainbow(nrow(peaks), start=0, end=4/6))[rank(peaks$pvalue_MAnorm)], pch='.', cex=.5, xlab="A", ylab="M", xlim=range(peaks$M), main=paste("H3K4me3 peaks btw HD and Control (n=",nrow(peaks),")",sep=""))
image(1,1:nrow(peaks),matrix(data=1:nrow(peaks),ncol=nrow(peaks),nrow=1),col=rev(rainbow(nrow(peaks), start=0, end=4/6)), xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n", mar=c(4,0,2,2))
axis(2,seq(1, sum(peaks$pvalue_MAnorm!=Inf & !is.na(peaks$pvalue_MAnorm)),length.out=length(seq(0,max(peaks$pvalue_MAnorm[peaks$pvalue_MAnorm!=Inf], na.rm=T),20))),seq(0,max(peaks$pvalue_MAnorm[peaks$pvalue_MAnorm!=Inf], na.rm=T),20))
mtext('-log10(pvalue)', side=1, line=1, adj=1)

peaks1=subset(peaks, FDR<10 & fold_enrichment>4)
par(mar=c(4,4,4,1))
layout(matrix(seq(2),nrow=1,ncol=2),widths=c(3,0.5),heights=1)
plot(peaks1$A_rescaled, peaks1$M_rescaled, col=rev(rainbow(nrow(peaks1), start=0, end=4/6))[rank(peaks1$pvalue_MAnorm)], pch='.', cex=.5, xlab="A", ylab="M", xlim=range(peaks$M), main=paste("Significant H3K4me3 peaks btw HD and Control (n=",nrow(peaks1),")\nSignificant peaks: MACS_FDR<10% and fold_enrichment >4",sep=""))
image(1,1:nrow(peaks1),matrix(data=1:nrow(peaks1),ncol=nrow(peaks1),nrow=1),col=rev(rainbow(nrow(peaks1), start=0, end=4/6)), xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n", mar=c(4,0,2,2))
axis(2,seq(1, sum(peaks1$pvalue_MAnorm!=Inf & !is.na(peaks1$pvalue_MAnorm)),length.out=length(seq(0,max(peaks1$pvalue_MAnorm[peaks1$pvalue_MAnorm!=Inf], na.rm=T),20))),seq(0,max(peaks1$pvalue_MAnorm[peaks1$pvalue_MAnorm!=Inf], na.rm=T),20))
mtext('-log10(pvalue)', side=1, line=1, adj=1)

# proxiaml vs. distal
par(mar=c(4,4,4,1))
layout(matrix(seq(2),nrow=1,ncol=2),widths=c(3,0.5),heights=1)
plot(peaks1$A_rescaled, peaks1$M_rescaled, col=ifelse(peaks1$inpromoter=="proximal", '#00000088','#ff000088'), pch='.', cex=.5, xlab="A", ylab="M", xlim=range(peaks$M), main=paste("Significant H3K4me3 peaks btw HD and Control (n=",nrow(peaks1),")\nSignificant peaks: MACS_FDR<10% and fold_enrichment >4",sep=""))
n=c(sum(peaks1$inpromoter=="proximal"), sum(peaks1$inpromoter!="proximal"))
legend("topleft", paste(c("proximal peaks","distal peaks"), " (n=", n,")", sep=""), col=c('#00000088','#ff000088'), pch=20, bty='n')
image(1,1:nrow(peaks1),matrix(data=1:nrow(peaks1),ncol=nrow(peaks1),nrow=1),col=rev(rainbow(nrow(peaks1), start=0, end=4/6)), xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n", mar=c(4,0,2,2))
axis(2,seq(1, sum(peaks1$pvalue_MAnorm!=Inf & !is.na(peaks1$pvalue_MAnorm)),length.out=length(seq(0,max(peaks1$pvalue_MAnorm[peaks1$pvalue_MAnorm!=Inf], na.rm=T),20))),seq(0,max(peaks1$pvalue_MAnorm[peaks1$pvalue_MAnorm!=Inf], na.rm=T),20))
mtext('-log10(pvalue)', side=1, line=1, adj=1)

# proximal peaks into different peak type (unique vs. common)
p1=subset(peaks1, inpromoter=="proximal")
par(mar=c(4,4,4,1))
layout(matrix(seq(2),nrow=1,ncol=2),widths=c(3,0.5),heights=1)
plot(subset(p1, peak_type=="common_peak1", c(A_rescaled, M_rescaled)), col='#0000aa66', pch=20, cex=.5, xlab="A", ylab="M", xlim=range(peaks$M), main=paste("Proximal H3K4me3 peaks (n=",nrow(p1),")\nSignificant peaks: MACS_FDR<10% and fold_enrichment >4",sep=""))
points(subset(p1, peak_type=="common_peak2", c(A_rescaled, M_rescaled)), col='#0000ff66', pch=20, cex=.5)
points(subset(p1, peak_type=="unique_peak1", c(A_rescaled, M_rescaled)), col='#ff000088', pch=20, cex=.5)
points(subset(p1, peak_type=="unique_peak2", c(A_rescaled, M_rescaled)), col='#00ff0088', pch=20, cex=.5)
n=table(p1$peak_type)
legend("topleft", paste(names(n), " (n=", n,")", sep=""), col=c('#0000aa66','#0000ff66', '#ff000088','#00ff0088'), pch=20, bty='n')
image(1,1:nrow(peaks1),matrix(data=1:nrow(peaks1),ncol=nrow(peaks1),nrow=1),col=rev(rainbow(nrow(peaks1), start=0, end=4/6)), xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n", mar=c(4,0,2,2))
axis(2,seq(1, sum(peaks1$pvalue_MAnorm!=Inf & !is.na(peaks1$pvalue_MAnorm)),length.out=length(seq(0,max(peaks1$pvalue_MAnorm[peaks1$pvalue_MAnorm!=Inf], na.rm=T),20))),seq(0,max(peaks1$pvalue_MAnorm[peaks1$pvalue_MAnorm!=Inf], na.rm=T),20))
mtext('-log10(pvalue)', side=1, line=1, adj=1)

# distal peaks into different peak type (unique vs. common)
p2=subset(peaks1, inpromoter=="distal")
par(mar=c(4,4,4,1))
layout(matrix(seq(2),nrow=1,ncol=2),widths=c(3,0.5),heights=1)
plot(subset(p2, peak_type=="common_peak1", c(A_rescaled, M_rescaled)), col='#0000aa66', pch=20, cex=.5, xlab="A", ylab="M", xlim=range(peaks$M), main=paste("Distal H3K4me3 peaks (n=",nrow(p2),")\nSignificant peaks: MACS_FDR<10% and fold_enrichment >4",sep=""))
points(subset(p2, peak_type=="common_peak2", c(A_rescaled, M_rescaled)), col='#0000ff66', pch=20, cex=.5)
points(subset(p2, peak_type=="unique_peak1", c(A_rescaled, M_rescaled)), col='#ff000088', pch=20, cex=.5)
points(subset(p2, peak_type=="unique_peak2", c(A_rescaled, M_rescaled)), col='#00ff0088', pch=20, cex=.5)
n=table(p2$peak_type)
legend("topleft", paste(names(n), " (n=", n,")", sep=""), col=c('#0000aa66','#0000ff66', '#ff000088','#00ff0088'), pch=20, bty='n')
image(1,1:nrow(peaks1),matrix(data=1:nrow(peaks1),ncol=nrow(peaks1),nrow=1),col=rev(rainbow(nrow(peaks1), start=0, end=4/6)), xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n", mar=c(4,0,2,2))
axis(2,seq(1, sum(peaks1$pvalue_MAnorm!=Inf & !is.na(peaks1$pvalue_MAnorm)),length.out=length(seq(0,max(peaks1$pvalue_MAnorm[peaks1$pvalue_MAnorm!=Inf], na.rm=T),20))),seq(0,max(peaks1$pvalue_MAnorm[peaks1$pvalue_MAnorm!=Inf], na.rm=T),20))
mtext('-log10(pvalue)', side=1, line=1, adj=1)


# proximal peaks into different peak type (unique vs. common)
p1=subset(peaks1, inpromoter=="proximal")
par(mar=c(4,4,4,1))
layout(matrix(seq(2),nrow=1,ncol=2),widths=c(3,0.5),heights=1)
plot(subset(p1, peak_type=="common_peak1", c(A_rescaled, M_rescaled)), col='#0000aa66', pch=20, cex=(p1$summit_height[p1$peak_type=="common_peak1"]/200)^2+0.3, xlab="A", ylab="M", xlim=range(peaks$M), main=paste("Proximal H3K4me3 peaks (n=",nrow(p1),")\nSignificant peaks: MACS_FDR<10% and fold_enrichment >4",sep=""))
points(subset(p1, peak_type=="common_peak2", c(A_rescaled, M_rescaled)), col='#0000ff66', pch=20, cex=(p1$summit_height[p1$peak_type=="common_peak2"]/200)^2+0.3)
points(subset(p1, peak_type=="unique_peak1", c(A_rescaled, M_rescaled)), col='#ff000088', pch=20, cex=(p1$summit_height[p1$peak_type=="unique_peak1"]/200)^2+0.3)
points(subset(p1, peak_type=="unique_peak2", c(A_rescaled, M_rescaled)), col='#00ff0088', pch=20, cex=(p1$summit_height[p1$peak_type=="unique_peak2"]/200)^2+0.3)
n=table(p1$peak_type)
legend("topleft", paste(names(n), " (n=", n,")", sep=""), col=c('#0000aa66','#0000ff66', '#ff000088','#00ff0088'), pch=20, bty='n')
legend("bottomright", title="H3K4me3 peak height", title.adj=0.5, legend=c(50,100,200,400), pt.cex=(c(50,100,200,400)/200)^2+0.3, col="#999999", pch=20,  horiz=T, bty='n')
image(1,1:nrow(peaks1),matrix(data=1:nrow(peaks1),ncol=nrow(peaks1),nrow=1),col=rev(rainbow(nrow(peaks1), start=0, end=4/6)), xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n", mar=c(4,0,2,2))
axis(2,seq(1, sum(peaks1$pvalue_MAnorm!=Inf & !is.na(peaks1$pvalue_MAnorm)),length.out=length(seq(0,max(peaks1$pvalue_MAnorm[peaks1$pvalue_MAnorm!=Inf], na.rm=T),20))),seq(0,max(peaks1$pvalue_MAnorm[peaks1$pvalue_MAnorm!=Inf], na.rm=T),20))
mtext('-log10(pvalue)', side=1, line=1, adj=1)

# distal peaks into different peak type (unique vs. common)
p2=subset(peaks1, inpromoter=="distal")
par(mar=c(4,4,4,1))
layout(matrix(seq(2),nrow=1,ncol=2),widths=c(3,0.5),heights=1)
plot(subset(p2, peak_type=="common_peak1", c(A_rescaled, M_rescaled)), col='#0000aa66', pch=20, cex=(p2$summit_height[p2$peak_type=="common_peak1"]/200)^2+0.3, xlab="A", ylab="M", xlim=range(peaks$M), main=paste("Distal H3K4me3 peaks (n=",nrow(p2),")\nSignificant peaks: MACS_FDR<10% and fold_enrichment >4",sep=""))
points(subset(p2, peak_type=="common_peak2", c(A_rescaled, M_rescaled)), col='#0000ff66', pch=20, cex=(p2$summit_height[p2$peak_type=="common_peak2"]/200)^2+0.3)
points(subset(p2, peak_type=="unique_peak1", c(A_rescaled, M_rescaled)), col='#ff000088', pch=20, cex=(p2$summit_height[p2$peak_type=="unique_peak1"]/200)^2+0.3)
points(subset(p2, peak_type=="unique_peak2", c(A_rescaled, M_rescaled)), col='#00ff0088', pch=20, cex=(p2$summit_height[p2$peak_type=="unique_peak2"]/200)^2+0.3)
n=table(p2$peak_type)
legend("topleft", paste(names(n), " (n=", n,")", sep=""), col=c('#0000aa66','#0000ff66', '#ff000088','#00ff0088'), pch=20, bty='n')
legend("bottomright", title="H3K4me3 peak height", title.adj=0.5, legend=c(50,100,200,400), pt.cex=(c(50,100,200,400)/200)^2+0.3, col="#999999", pch=20, horiz=T, bty='n')
image(1,1:nrow(peaks1),matrix(data=1:nrow(peaks1),ncol=nrow(peaks1),nrow=1),col=rev(rainbow(nrow(peaks1), start=0, end=4/6)), xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n", mar=c(4,0,2,2))
axis(2,seq(1, sum(peaks1$pvalue_MAnorm!=Inf & !is.na(peaks1$pvalue_MAnorm)),length.out=length(seq(0,max(peaks1$pvalue_MAnorm[peaks1$pvalue_MAnorm!=Inf], na.rm=T),20))),seq(0,max(peaks1$pvalue_MAnorm[peaks1$pvalue_MAnorm!=Inf], na.rm=T),20))
mtext('-log10(pvalue)', side=1, line=1, adj=1)

dev.off()
