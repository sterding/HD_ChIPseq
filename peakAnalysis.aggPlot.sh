## script to generate data to draw aggregation plot for enhancer regions 

ANNOTATION=$GENOME/Annotation/Genes

cd ~/projects/projects/HD/data

#============================================================
# Enhancers feature data. 
#============================================================
# CAGE
wget http://fantom.gsc.riken.jp/5/datahub/hg19/reads/ctssTotalCounts.fwd.bw
wget http://fantom.gsc.riken.jp/5/datahub/hg19/reads/ctssTotalCounts.rev.bw

ln -s ctssTotalCounts.fwd.bw CAGE.fwd.bigwig
ln -s ctssTotalCounts.rev.bw CAGE.rev.bigwig

# histone marks
#http://www.ncbi.nlm.nih.gov/geo/roadmap/epigenomics/?view=samples&sample=brain%2C%20mid%20frontal%20lobe
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM773nnn/GSM773014/suppl/GSM773014_BI.Brain_Mid_Frontal_Lobe.H3K4me1.149.wig.gz
wget http://www.genboree.org/EdaccData/Current-Release/sample-experiment/Brain_Mid_Frontal_Lobe/Histone_H3K4me1/BI.Brain_Mid_Frontal_Lobe.H3K4me1.112.wig.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM773nnn/GSM773012/suppl/GSM773012_BI.Brain_Mid_Frontal_Lobe.H3K4me3.149.wig.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1112nnn/GSM1112810/suppl/GSM1112810_BI.Brain_Mid_Frontal_Lobe.H3K27ac.112.wig.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM772nnn/GSM772833/suppl/GSM772833_BI.Brain_Mid_Frontal_Lobe.H3K27me3.149.wig.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM773nnn/GSM773013/suppl/GSM773013_BI.Brain_Mid_Frontal_Lobe.H3K36me3.149.wig.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM670nnn/GSM670021/suppl/GSM670021_BI.Brain_Mid_Frontal_Lobe.H3K9ac.112.wig.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM772nnn/GSM772834/suppl/GSM772834_BI.Brain_Mid_Frontal_Lobe.H3K9me3.149.wig.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM916nnn/GSM916085/suppl/GSM916085_BI.Brain_Mid_Frontal_Lobe.RRBS.149.wig.gz

gunzip *.wig.gz

for i in H3K4me1 H3K4me3 H3K27ac; do wigToBigWig *$i*.wig $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size $i.bigwig & done

ln -s BI.Brain_Mid_Frontal_Lobe.H3K4me1.112.bigWig H3K4me1.bigwig 

#TFBS
ln -s ~/neurogen/rnaseq_PD/results/eRNA/externalData/TFBS/wgEncodeRegTfbsClusteredV3.bw TFBS.bigwig
#DNase
wget https://www.encodedcc.org/files/ENCFF000SFV/@@download/ENCFF000SFV.bigWig
ln -s ~/neurogen/rnaseq_PD/results/eRNA/externalData/DNase/UW.Fetal_Brain.ChromatinAccessibility.H-24510.DNase.DS20780.bw DNase.bigwig

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeOpenChromDnase/wgEncodeOpenChromDnaseFrontalcortexocPk.narrowPeak.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeOpenChromDnase/wgEncodeOpenChromDnaseFrontalcortexocSig.bigWig
#wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeOpenChromDnase/wgEncodeOpenChromDnaseFrontalcortexocBaseOverlapSignal.bigWig

#Conservation
ln -s ~/neurogen/rnaseq_PD/results/eRNA/externalData/Conservation/vertebrate/phyloP46way.wigFix.bigwig Conservation.bigwig

#============================================================
# get bin signal
#============================================================
mkdir peaks.all peaks.proximal peaks.distal
#extend to 1000bp from the summit of DNase narrowpeak (if any) or the center points of a H3k4me3 peak
more +2 h3k4me3_peaks_union.signal.neighborhood.annotated.bed.v3.xls | cut -f2- | intersectBed -a stdin -b wgEncodeOpenChromDnaseFrontalcortexocPk.narrowPeak -wao | sort -k1,1 -k2,2n -k21,21gr  | awk '{OFS="\t"; if($2!=id) {print; id=$2}}' | awk '{OFS="\t"; mid=($16!=-1)?($16+$24):int(($2+$3)/2); print $1,mid-500,mid+500;}' > peaks.all/h3k4me3_peaks_union.1000bp.bed
grep -w distal h3k4me3_peaks_union.signal.neighborhood.annotated.bed.v3.xls | cut -f2- | intersectBed -a stdin -b wgEncodeOpenChromDnaseFrontalcortexocPk.narrowPeak -wao | sort -k1,1 -k2,2n -k21,21gr  | awk '{OFS="\t"; if($2!=id) {print; id=$2}}' | awk '{OFS="\t"; mid=($16!=-1)?($16+$24):int(($2+$3)/2); print $1,mid-500,mid+500;}' > peaks.distal/h3k4me3_peaks_union.1000bp.bed
grep -w proximal h3k4me3_peaks_union.signal.neighborhood.annotated.bed.v3.xls | cut -f2- | intersectBed -a stdin -b wgEncodeOpenChromDnaseFrontalcortexocPk.narrowPeak -wao | sort -k1,1 -k2,2n -k21,21gr  | awk '{OFS="\t"; if($2!=id) {print; id=$2}}' | awk '{OFS="\t"; mid=($16!=-1)?($16+$24):int(($2+$3)/2); print $1,mid-500,mid+500;}' > peaks.proximal/h3k4me3_peaks_union.1000bp.bed

#
#more +2 h3k4me3_peaks_union.signal.neighborhood.annotated.bed.v3.xls | awk '{OFS="\t"; mid=int(($3+$4)/2); print $2,mid-500,mid+500;}' > peaks.all/h3k4me3_peaks_union.bed
#grep -w distal h3k4me3_peaks_union.signal.neighborhood.annotated.bed.v3.xls | awk '{OFS="\t"; mid=int(($3+$4)/2); print $2,mid-500,mid+500;}' > peaks.distal/h3k4me3_peaks_union.bed
#grep -w proximal h3k4me3_peaks_union.signal.neighborhood.annotated.bed.v3.xls | awk '{OFS="\t"; mid=int(($3+$4)/2); print $2,mid-500,mid+500;}' >peaks.proximal/h3k4me3_peaks_union.bed
ls *.bigwig | parallel '~/neurogen/pipeline/RNAseq/bin/toBinRegionsOnBigwig.sh {} peaks.all/h3k4me3_peaks_union.1000bp.bed 200 > peaks.all/{}.bins' &
ls *.bigwig | parallel '~/neurogen/pipeline/RNAseq/bin/toBinRegionsOnBigwig.sh {} peaks.distal/h3k4me3_peaks_union.1000bp.bed 200 > peaks.distal/{}.bins' &
ls *.bigwig | parallel '~/neurogen/pipeline/RNAseq/bin/toBinRegionsOnBigwig.sh {} peaks.proximal/h3k4me3_peaks_union.1000bp.bed 200 > peaks.proximal/{}.bins' &

#============================================================
# draw aggregation plot (to run in R console)
#============================================================
Rscript ~/projects/HD/src/peakAnalysis.aggPlot.R ~/projects/HD/data/peaks.all ~/projects/HD/aggregation.peaks.all.pdf
Rscript ~/projects/HD/src/peakAnalysis.aggPlot.R ~/projects/HD/data/peaks.proximal ~/projects/HD/aggregation.peaks.proximal.pdf
Rscript ~/projects/HD/src/peakAnalysis.aggPlot.R ~/projects/HD/data/peaks.distal ~/projects/HD/aggregation.peaks.distal.pdf
