#######################################
# Script to call gene/tx with differential H3k4me3 signal at promoter regions, using DESeq2
# Author: Xianjun Dong
# Date: 2013-12-14
# Version: 0.0
#######################################

#!/bin/sh


# required files
annotation=$GENOME/hg19/Annotation/Genes/gencode.v17.annotation.gtf.gz
promoters=$annotation.promoter_4k_round_TSS.bed
# #zcat $annotation | fgrep -w transcript | sed 's/[";]//g;' | awk '{OFS="\t"; tss=($7=="+")?($4-1):($5-1); if(tss<2000) tss=2000; print $1, tss-2000,tss+2000,$12"___"$10"___"$14"___"$18, 4000, $7}'  > $annotation.promoter_4k_round_TSS.bed

echo "Step 1: calculate reads count of h3k4me3 in the promoter regions"
#######################################################################

for i in ~/scratch/combined_Lane*/accepted_hits.bam ~/scratch/h3k4me3_HD_batch1*/accepted_hits.bam ~/scratch/Schraham_*/accepted_hits.bam; do
    qsub ~/projects/bu_neuro/src/DE_ChIPseq.sh $i $promoters
done

## wait until all jobs are done

echo "Step 2: cat the readscount file together"
#######################################################################

RC_in_promoter=~/projects/bu_neuro/data/h3k4me3.in.promoter.allsamples.txt
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










