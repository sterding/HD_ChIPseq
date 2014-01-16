# Run this R script after getting reads count for all features in all samples  [deg_sge.sh]
# paste ~/scratch/HD_*/*.tab | cut -f1,2,4,6,8,10,12,14 > ../results/tophat.mm2.accepted_hits.hg19.sam.hgseqcount.tab

library(DESeq)
countsTable <- read.delim( "../results/DESeq/tophat.mm2.accepted_hits.hg19.sam.hgseqcount.tab", header=T, stringsAsFactors=TRUE )
geneanno = read.delim("../data/Ensembl59.transcripts.annotation.tab")  # file download from Biomart for Ens v59

rownames( countsTable ) <- countsTable$ID
countsTable <- countsTable[ , -1 ]
conds <- substring(colnames(countsTable), 1, 2)

cds <- newCountDataSet( countsTable[grep("ENST", rownames(countsTable)),], conds )

libsizes <- apply(countsTable, 2, sum)
sizeFactors(cds) <- libsizes
cds <- estimateSizeFactors( cds )
cds <- estimateVarianceFunctions( cds )
res <- nbinomTest( cds, "Ct", "HD")

sig=res[res$padj<0.1 & !is.na(res$padj),]

# top ones
sigsort=sig[with(sig, order(padj)), ]

# which genes are they?
sigsortgenes = cbind(sigsort[match(intersect(sigsort$id, geneanno$Trans.ID), sigsort$id),],
                     geneanno[match(intersect(sigsort$id, geneanno$Trans.ID), geneanno$Trans.ID),c(1,2,8:11)])

write(head(sigsortgenes, n=200)$Gene.ID, "../results/DESeq/GENEID.TOP.txt")


## -------------  deal with cuffdiff's output -------------
fpkm0 = read.delim( "../results/HD_cuffdiff/genes.fpkm_tracking", header=T, stringsAsFactors=TRUE )
fpkm = fpkm0[,grep("gene|FPKM", colnames(fpkm0))]
rownames( fpkm ) <- fpkm0$tracking_id
fpkm <- fpkm[ , -1 ]

# remove all-0 records
fpkm=fpkm[apply(fpkm, 1, sum)>0,]  # 12193/48065 = 25% left
fpkm[fpkm==0]=0.001

library(gplots)
#pdf(file="RNAseq.K562.RPKM.heatmap.2.pdf", height =10, width = 8)
heatmap.2(as.matrix(log(fpkm[1:1000,])),
          col=greenred, #colorpanel(8,"green","white","red"),
          trace='none', margins=c(8, 8),
          key=T, keysize = 1, scale='none', cexRow=0.5,cexCol=0.6, density.info= 'none',
          main = "RNAseq.FPKM.HD-Control.heatmap",
           xlab = "Experiments (HD: Huntington's Disease; Ct: Control)",
           ylab = "Transcript ID")
