## script to draw aggregation plot (to run in R console)
## usage: Rscript ~/projects/HD/src/peakAnalysis.aggPlot.R ~/projects/HD/data/peaks.all ~/projects/HD/aggregation.peaks.pdf
#============================================================
args<-commandArgs(TRUE)
input_dir=args[1] #"~/projects/HD/data/peaks.all";
output_pdf=args[2] # "~/projects/HD/aggregation.peaks.pdf"

pdf(output_pdf, width=3, height=8)
par(mfcol=c(5,1), mar = c(2,4,1,1))

#CAGE
df1=read.table(paste(input_dir, "CAGE.fwd.bigwig.bins", sep="/"), header=F)[,-1]
df2=read.table(paste(input_dir, "CAGE.rev.bigwig.bins", sep="/"), header=F)[,-1]
df=cbind(df1, df2);
df=df[apply(df,1,sum)>0,];
m=nrow(df);
df=apply(df, 2, function(x) mean(x, trim=0.001)); N=length(df)/2
plot(df[1:N], type="l", col="green", ylab=NA, xlab=NA, xaxt="n", main="", ylim=range(df)) 
points(df[-c(1:N)], type="l", col="blue")
legend("topright", c("CAGE +","CAGE -"), col=c("green","blue"), lty=1, bty='n')
axis(1, at=50*c(0:4), labels=c("-500 bp","","0","","500 bp"))
mtext(side = 2, "Mean signal of CAGE\n (total ctss counts in FANTOM5)", line = 2, cex=.6)

#Histone
df1=read.table(paste(input_dir, "H3K27ac.bigwig.bins", sep="/"), header=F)[,-1]
df2=read.table(paste(input_dir, "H3K4me1.bigwig.bins", sep="/"), header=F)[,-1]
df3=read.table(paste(input_dir, "H3K4me3.bigwig.bins", sep="/"), header=F)[,-1]
df4=read.table(paste(input_dir, "H3K9ac.bigwig.bins", sep="/"), header=F)[,-1]

df=cbind(df1, df2, df3, df4)
df=df[apply(df,1,sum)>0,];
m=nrow(df);
df=apply(df, 2, function(x) mean(x, trim=0.01)); N=length(df)/4
plot(df[1:N], type="l", col="green", ylab=NA, xlab=NA, xaxt="n", ylim=c(0,40)) #range(df)) 
points(df[N+c(1:N)], type="l", col="blue")
points(df[2*N+c(1:N)], type="l", col="red")
points(df[3*N+c(1:N)], type="l", col="gold")
legend("topright", c("H3K4me3","H3K27ac","H3K4me1", "H3K9ac"), col=c("red","green","blue","gold"), lty=1, bty='n');
mtext(side = 2, "Mean signal of histone marks \n in Brodmann area 9/46", line = 2, cex=.6)
axis(1, at=50*c(0:4), labels=c("-500 bp","","0","","500 bp"))

#TFBS
df=read.table(paste(input_dir, "TFBS.bigwig.bins", sep="/"), header=F)[,-1]
df=df[apply(df,1,sum)>0,];
m=nrow(df);
df=apply(df, 2, function(x) mean(x, trim=0.01));
plot(df, type="l", col="blue", ylab=NA, xlab=NA, xaxt="n") 
axis(1, at=50*c(0:4), labels=c("-500 bp","","0","","500 bp"))
mtext(side = 2, "Mean TF binding occurances \n of 161 ENCODE TFs", line = 2, cex=.6)

#Conservation
df=read.table(paste(input_dir, "Conservation.bigwig.bins", sep="/"), header=F)[,-1]
df=df[apply(df,1,sum)>0,];
m=nrow(df);
df=apply(df, 2, function(x) mean(x, trim=0.01));
plot(df, type="l", col="blue", ylab=NA, xlab=NA, xaxt="n", ylim=range(df,0.2)) 
axis(1, at=50*c(0:4), labels=c("-500 bp","","0","","500 bp"))
mtext(side = 2, "Mean phyloP score \nin 46-way conservation", line = 2, cex=.6)

#DNase
df=read.table(paste(input_dir, "DNase.bigwig.bins", sep="/"), header=F)[,-1]
df=df[apply(df,1,sum)>0,];
m=nrow(df);
df=apply(df, 2, function(x) mean(x, trim=0.01));
plot(df, type="l", col="blue", ylab=NA, xlab=NA, xaxt="n") 
axis(1, at=50*c(0:4), labels=c("-500 bp","","0","","500 bp"))
mtext(side = 2, "Mean DNase signal \nin fetal brain (Roadmap Epigenomics)", line = 2, cex=.6)
        
dev.off()