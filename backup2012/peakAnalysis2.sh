#!/bin/bash

# required files
annotation=$GENOME/hg19/Annotation/Genes/gencode.v17.annotation.gtf.gz
#zcat $annotation | fgrep -w transcript | sed 's/[";]//g;' | awk '{OFS="\t"; tss=($7=="+")?($4-1):($5-1); if(tss<2000) tss=2000; print $1, tss-2000,tss+2000,$12"___"$10"___"$14"___"$18, 4000, $7}'  > $annotation.promoter_4k_round_TSS.bed

#expression=$HOME/projects/bu_neuro/data/mlhd_DESeq_normalized_counts_all_controls.txt  # V1 (21HD and 54Ct)
expression=$HOME/projects/bu_neuro/data/mlhd_DESeq_normalized_counts_all_controls_strict.txt  # V2 (21HD and 50Ct, and using strict DEseq criteria)

# input
MAnorm_result=$1  #MAnorm_result.xls

echo "Step 1: define peaks into proximal/distal"

echo `head -n1 ../HD.vs.Ct_M_g43yr/MAnorm_result.xls` inpromoter | sed 's/ /\t/g;s/#//g' > $MAnorm_result.inpromoter.tab
bedtools intersect -b $annotation.promoter_4k_round_TSS.bed -a <(grep -P "^chr[^\s]" $MAnorm_result) -wa -u | awk '{OFS="\t"; print $0, "proximal"}' >> $MAnorm_result.inpromoter.tab
bedtools intersect -b $annotation.promoter_4k_round_TSS.bed -a <(grep -P "^chr[^\s]" $MAnorm_result) -wa -v | awk '{OFS="\t"; print $0, "distal"}' >> $MAnorm_result.inpromoter.tab

echo "Step 2: genes with H3k4kem3 peaks in proximal promoter"

## assign the closest peak to each TSS first, and then assigning the overall gene expression to the strongest peak in ALL isoforms' promoter.
bedtools intersect -a $annotation.promoter_4k_round_TSS.bed -b <(grep -P "^chr[^\s]" $MAnorm_result) -wao | awk '{OFS="\t"; dis_tss2sumit=($8+$9)/2-($2+$3)/2; if(dis_tss2sumit<0) dis_tss2sumit=-dis_tss2sumit; print $0, dis_tss2sumit}' | sort -k4,4 -k17,17n | awk '{if($4!=id) {id=$4; print;}}' | sed 's/___/\t/g' | sort -k5,5 -k17,17nr | awk '{if($5!=id) {id=$5; print;}}' > allgenes.H3k4me3inpromoter.tab

echo "Step 3: merge with expression level"
# merge with expression level

join -1 5 -2 1 -a 1 allgenes.H3k4me3inpromoter.tab <(sed 's/"//g' $expression | sort -k1,1) | sed 's/ /\t/g' > allgenes.H3k4me3inpromoter.expression.tab

intersectBed -a /home/dongx/nearline/genomes/hg19/Annotation/Genes/gencode.v17.annotation.gtf.gz.promoter_4k_round_TSS.bed -b <(fgrep -w proximal MAnorm_result.xls.inpromoter.tab | grep common) -u | sed 's/___/\t/g' | cut -f5 | sort -u > id.common.txt
fgrep -f id.common.txt $expression | 