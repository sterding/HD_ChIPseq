args<-commandArgs(TRUE)
peakfile=args[1]  # MAnorm_result.xls + column of promoter/enhancer info

peaks=read.table(peakfile, header=T)
#columns: chr	start	end	description	raw_read_1	raw_read_2	M_value_rescaled	A_value_rescaled	-log10(p-value)	inpromoter
names(peaks)[9]="pvalue_MAnorm"

pdf("peakAnalysis.pdf", width=8, height=6, colormodel='cmyk')
par(mar=c(4,4,4,1))
layout(matrix(seq(2),nrow=1,ncol=2),widths=c(3,0.5),heights=1)
plot(peaks$A_value_rescaled, peaks$M_value_rescaled, col=rev(rainbow(nrow(peaks), start=0, end=4/6))[rank(peaks$pvalue_MAnorm)], pch='.', cex=.5, xlab="A", ylab="M", xlim=range(peaks$M), main=paste("H3K4me3 peaks btw HD and Control (n=",nrow(peaks),")",sep=""))
image(1,1:nrow(peaks),matrix(data=1:nrow(peaks),ncol=nrow(peaks),nrow=1),col=rev(rainbow(nrow(peaks), start=0, end=4/6)), xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n", mar=c(4,0,2,2))
axis(2,seq(1, sum(peaks$pvalue_MAnorm!=Inf & !is.na(peaks$pvalue_MAnorm)),length.out=length(seq(0,max(peaks$pvalue_MAnorm[peaks$pvalue_MAnorm!=Inf], na.rm=T),20))),seq(0,max(peaks$pvalue_MAnorm[peaks$pvalue_MAnorm!=Inf], na.rm=T),20))
mtext('-log10(pvalue)', side=1, line=1, adj=1)

peaks1=peaks

# proxiaml vs. distal
par(mar=c(4,4,4,1))
layout(matrix(seq(2),nrow=1,ncol=2),widths=c(3,0.5),heights=1)
plot(peaks1$A_value_rescaled, peaks1$M_value_rescaled, col=ifelse(peaks1$inpromoter=="proximal", '#00000088','#ff000088'), pch='.', cex=.5, xlab="A", ylab="M", xlim=range(peaks$M), main=paste("Significant H3K4me3 peaks btw HD and Control (n=",nrow(peaks1),")\nSignificant peaks: MACS_FDR<10% and fold_enrichment >4",sep=""))
n=c(sum(peaks1$inpromoter=="proximal"), sum(peaks1$inpromoter!="proximal"))
legend("topleft", paste(c("proximal peaks","distal peaks"), " (n=", n,")", sep=""), col=c('#00000088','#ff000088'), pch=20, bty='n')
image(1,1:nrow(peaks1),matrix(data=1:nrow(peaks1),ncol=nrow(peaks1),nrow=1),col=rev(rainbow(nrow(peaks1), start=0, end=4/6)), xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n", mar=c(4,0,2,2))
axis(2,seq(1, sum(peaks1$pvalue_MAnorm!=Inf & !is.na(peaks1$pvalue_MAnorm)),length.out=length(seq(0,max(peaks1$pvalue_MAnorm[peaks1$pvalue_MAnorm!=Inf], na.rm=T),20))),seq(0,max(peaks1$pvalue_MAnorm[peaks1$pvalue_MAnorm!=Inf], na.rm=T),20))
mtext('-log10(pvalue)', side=1, line=1, adj=1)

# proximal peaks into different peak type (unique vs. common)
p1=subset(peaks1, inpromoter=="proximal")
par(mar=c(4,4,4,1))
layout(matrix(seq(2),nrow=1,ncol=2),widths=c(3,0.5),heights=1)
plot(subset(p1, description=="common_peak1", c(A_value_rescaled, M_value_rescaled)), col='#0000aa66', pch=20, cex=.5, xlab="A", ylab="M", xlim=range(peaks$M), main=paste("Proximal H3K4me3 peaks (n=",nrow(p1),")\nSignificant peaks: MACS_FDR<10% and fold_enrichment >4",sep=""))
points(subset(p1, description=="common_peak2", c(A_value_rescaled, M_value_rescaled)), col='#0000ff66', pch=20, cex=.5)
points(subset(p1, description=="unique_peak1", c(A_value_rescaled, M_value_rescaled)), col='#ff000088', pch=20, cex=.5)
points(subset(p1, description=="unique_peak2", c(A_value_rescaled, M_value_rescaled)), col='#00ff0088', pch=20, cex=.5)
n=table(p1$description)
legend("topleft", paste(names(n), " (n=", n,")", sep=""), col=c('#0000aa66','#0000ff66', '#ff000088','#00ff0088'), pch=20, bty='n')
image(1,1:nrow(peaks1),matrix(data=1:nrow(peaks1),ncol=nrow(peaks1),nrow=1),col=rev(rainbow(nrow(peaks1), start=0, end=4/6)), xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n", mar=c(4,0,2,2))
axis(2,seq(1, sum(peaks1$pvalue_MAnorm!=Inf & !is.na(peaks1$pvalue_MAnorm)),length.out=length(seq(0,max(peaks1$pvalue_MAnorm[peaks1$pvalue_MAnorm!=Inf], na.rm=T),20))),seq(0,max(peaks1$pvalue_MAnorm[peaks1$pvalue_MAnorm!=Inf], na.rm=T),20))
mtext('-log10(pvalue)', side=1, line=1, adj=1)

# distal peaks into different peak type (unique vs. common)
p2=subset(peaks1, inpromoter=="distal")
par(mar=c(4,4,4,1))
layout(matrix(seq(2),nrow=1,ncol=2),widths=c(3,0.5),heights=1)
plot(subset(p2, description=="common_peak1", c(A_value_rescaled, M_value_rescaled)), col='#0000aa66', pch=20, cex=.5, xlab="A", ylab="M", xlim=range(peaks$M), main=paste("Distal H3K4me3 peaks (n=",nrow(p2),")\nSignificant peaks: MACS_FDR<10% and fold_enrichment >4",sep=""))
points(subset(p2, description=="common_peak2", c(A_value_rescaled, M_value_rescaled)), col='#0000ff66', pch=20, cex=.5)
points(subset(p2, description=="unique_peak1", c(A_value_rescaled, M_value_rescaled)), col='#ff000088', pch=20, cex=.5)
points(subset(p2, description=="unique_peak2", c(A_value_rescaled, M_value_rescaled)), col='#00ff0088', pch=20, cex=.5)
n=table(p2$description)
legend("topleft", paste(names(n), " (n=", n,")", sep=""), col=c('#0000aa66','#0000ff66', '#ff000088','#00ff0088'), pch=20, bty='n')
image(1,1:nrow(peaks1),matrix(data=1:nrow(peaks1),ncol=nrow(peaks1),nrow=1),col=rev(rainbow(nrow(peaks1), start=0, end=4/6)), xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n", mar=c(4,0,2,2))
axis(2,seq(1, sum(peaks1$pvalue_MAnorm!=Inf & !is.na(peaks1$pvalue_MAnorm)),length.out=length(seq(0,max(peaks1$pvalue_MAnorm[peaks1$pvalue_MAnorm!=Inf], na.rm=T),20))),seq(0,max(peaks1$pvalue_MAnorm[peaks1$pvalue_MAnorm!=Inf], na.rm=T),20))
mtext('-log10(pvalue)', side=1, line=1, adj=1)

dev.off()
