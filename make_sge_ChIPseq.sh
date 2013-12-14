## Step1. bash script for running chipseq pipeline on cluster
# After this, go to peakCalling.sh (Step2)

#!/bin/sh
#$ -V
#$ -pe openmpi 12
#$ -cwd
#$ -o $HOME/sge_jobs_output/sge_job.$JOB_ID.out -j y
#$ -S /bin/bash
#$ -l mem_free=5G

### for i in ~/nearline/BU/h3k4me3/chipseq_h3k4me3_May2013/130420-376_and_130429-380_Richard_Myers/*; do j=`basename $i`; qsub ~/projects/bu_neuro/src/make_sge_ChIPseq.sh $i ${j/.*/}; done
#Your job 9376710 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376711 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376712 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376713 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376714 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376715 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376716 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376717 ("make_sge_ChIPseq.sh") has been submitted
### for i in ~/nearline/BU/h3k4me3/chipseq_JF/*.gz; do j=`basename $i`; qsub ~/projects/bu_neuro/src/make_sge_ChIPseq.sh $i ${j/.*/}; done
#Your job 9376718 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376719 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376720 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376721 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376722 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376723 ("make_sge_ChIPseq.sh") has been submitted
#
#
## for i in ~/nearline/BU/h3k4me3/control_from_Schraham/*.fq; do j=`basename $i`; qsub ~/projects/bu_neuro/src/make_sge_ChIPseq.sh $i Schraham_${j/.*/}; done
#Your job 9376724 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376725 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376726 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376727 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376728 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376729 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376730 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376731 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376732 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376733 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376734 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376735 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376736 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376737 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376738 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376739 ("make_sge_ChIPseq.sh") has been submitted
#
## for i in ~/nearline/BU/h3k4me3/h3k4me3_HD_batch1/*.fq; do j=`basename $i`; qsub ~/projects/bu_neuro/src/make_sge_ChIPseq.sh $i h3k4me3_HD_batch1_${j/.*/}; done
#Your job 9376740 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376741 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376742 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376743 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376744 ("make_sge_ChIPseq.sh") has been submitted
#Your job 9376745 ("make_sge_ChIPseq.sh") has been submitted


###########################################
############## 1. Configuring
###########################################

index=hg19
adaptor='yes'

readsfile=$1
samplename=$2

controlIP=~/nearline/BU/h3k4me3/2037in.dt1.bed

[ -d ~/scratch/$samplename ] || mkdir -p ~/scratch/$samplename
[ -d ~/scratch/ucsc ] || mkdir ~/scratch/ucsc

export BOWTIE_INDEXES=$GENOME/$index/Sequence/BowtieIndex/
export BOWTIE2_INDEXES=$GENOME/$index/Sequence/Bowtie2Index/
export ANNOTATION=$GENOME/$index/Annotation/Genes
export SEQUENCE=$GENOME/$index/Sequence/WholeGenomeFasta

phred=`getphred $readsfile`; scoreoption=""; # default
[ "$phred" == "Phred+64" ] && scoreoption="--phred64";
[ "$phred" == "Solexa+64" ] && scoreoption="--solexa-quals";
echo "scoreoption: $scoreoption for $readsfile"

############################################
################  2. quality filter: adaptor removal/clip
############################################
#
#if [[ "$adaptor" == "yes" ]]; then
#    far -s $readsfile -t $readsfile.clipped -f $far -as $adaptorsequence --cut-off 5 --min-overlap 10  --min-readlength 20 --trim-end $trimend --adaptive-overlap yes --nr-threads 8 --max-uncalled 30
#    mv $readsfile $readsfile.orgin
#    mv $readsfile.clipped $readsfile
#fi
#
##############################################
################# 3. QC
#############################################
#
#fastqc --outdir ~/scratch/$samplename $(echo $readsfile | sed 's/,/ /g')
#rm ~/scratch/$samplename/*fastqc.zip\
## QC
#mv ~/scratch/$samplename/*_fastqc ~/scratch/ucsc
#
#
##############################################
################## 4. mapping to the genome
##############################################
##
#### using bowtie2
### [reporting in default mode] look for multiple alignments, report the best one, similar to -M 1 in bowtie
#echo "bowtie2 -x genome_offrate3 $scoreoption -p 12 --mm -U <(zcat $readsfile) -S ~/scratch/$samplename/accepted_hits.sam 2> ~/scratch/$samplename/mapping.summary"
#bowtie2 -x genome_offrate3 $scoreoption -p 12 --mm -U <(zcat $readsfile) -S ~/scratch/$samplename/accepted_hits.sam 2> ~/scratch/$samplename/mapping.summary

### using bowtie
### look for multiple alignments, report the best one; allow 1 mismatch
#echo "bowtie genome_offrate3 -p 12 --mm -v 1 -M 1 --best --strata -q $readsfile -S ~/scratch/$samplename/accepted_hits.sam 2> ~/scratch/$samplename/mapping.summary"
#[[ $readsfile =~ "fq$|fastq$" ]] && bowtie genome_offrate3 -p 12 -v 1 -M 1 --best --strata -q $readsfile -S ~/scratch/$samplename/accepted_hits.sam 2> ~/scratch/$samplename/mapping.summary
#[[ $readsfile =~ "gz$" ]] && bowtie genome_offrate3 -p 12 -v 1 -M 1 --best --strata -q <(zcat $readsfile) -S ~/scratch/$samplename/accepted_hits.sam 2> ~/scratch/$samplename/mapping.summary
#
############################################
################ 5. post-processing, format converting
############################################
##
### sam - bam - sorted - index
##
cd ~/scratch/$samplename
#
#samtools view -Sbut $BOWTIE2_INDEXES/genome.fai accepted_hits.sam | samtools sort - accepted_hits.sorted
#mv accepted_hits.sorted.bam accepted_hits.bam
#samtools index accepted_hits.bam

# faster version
split -l 2000000 -d accepted_hits.sam accepted_hits_sam
> paraFile;
for i in accepted_hits_sam*; do echo "sam2bed -v bed12=F -v XSstrand=F $i | sort-bed - > $i.bed" >>  paraFile; done
ParaFly -c paraFile -CPU 12
sort -k1,1 -k2,2n -m accepted_hits_sam*.bed > accepted_hits.bed
rm accepted_hits_sam*.bed
# classical version
#sam2bed -v bed12=F -v XSstrand=F accepted_hits.sam > accepted_hits.bed

sort -k1,1 accepted_hits.bed | bedItemOverlapCount $index -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n > accepted_hits.bedGraph
bedGraphToBigWig accepted_hits.bedGraph $ANNOTATION/ChromInfo.txt accepted_hits.bw


###########################################
############### 6. peakcalling
###########################################
cd ~/scratch/$samplename
rm -rf ${samplename}_MACS_wiggle
macs14 -t accepted_hits.bed -c $controlIP -f BED -n $samplename -w -S > macs14.log

###########################################
############## 6. prepare for tracks files to display on UCSC
###########################################

# make index for the (sorted) BAM
#cp accepted_hits.bam ~/scratch/ucsc/$samplename.accepted_hits.bam
#cp accepted_hits.bam.bai ~/scratch/ucsc/$samplename.accepted_hits.bam.bai
#cp accepted_hits.bw ~/scratch/ucsc/$samplename.accepted_hits.bw

echo "DONE!"