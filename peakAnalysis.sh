#!/bin/bash

cd /Users/xdong/projects/HD/data/scratch
# required files
annotation=../gencode.v17.annotation.gtf.gz
#zcat $annotation | fgrep -w transcript | sed 's/[";]//g;' | awk '{OFS="\t"; tss=($7=="+")?($4-1):($5-1); if(tss<2000) tss=2000; print $1, tss-2000,tss+2000,$12"___"$10"___"$14"___"$18, 4000, $7}'  > $annotation.promoter_4k_round_TSS.bed
# bedtools slop -b 18000 -i $annotation.promoter_4k_round_TSS.bed -g hg19.genome > $annotation.promoter_2M_round_TSS.bed

#expression=$HOME/projects/bu_neuro/data/mlhd_DESeq_normalized_counts_all_controls.txt  # V1 (21HD and 54Ct)
#expression=$HOME/projects/bu_neuro/data/mlhd_DESeq_normalized_counts_all_controls_strict.txt  # V2 (21HD and 50Ct, and using strict DEseq criteria)

# input
peakfile1=h3k4me3_peaks_intersected.Ct.bed  # from merged peaks of samples of same group
peakfile2=h3k4me3_peaks_intersected.HD.bed

echo "Step 1: classify peaks into common/uniq"

intersectBed -a $peakfile1 -b $peakfile2 -wao | awk '{OFS="\t"; $4=($4==".")?"uniq":"common"; print $1, $2, $3, $4;}' > $peakfile1.commonuniq
intersectBed -a $peakfile2 -b $peakfile1 -wao | awk '{OFS="\t"; $4=($4==".")?"uniq":"common"; print $1, $2, $3, $4;}' > $peakfile2.commonuniq

for i in *commonuniq; do echo $i; cut -f4 $i | sort | uniq -c; done
#h3k4me3_peaks_intersected.Ct.bed.commonuniq
#  27011 common
#   6133 uniq
#h3k4me3_peaks_intersected.HD.bed.commonuniq
#  23367 common
#   3941 uniq

echo "Step 2: classify peaks into proximal/distal"

for i in *commonuniq; do
    echo $i; 
    bedtools intersect -a $i -b $annotation.promoter_4k_round_TSS.bed -wa -u | cut -f1-9,14,18-23 | awk '{OFS="\t"; print $0, "proximal"}' > $i.inpromoter
    bedtools intersect -a $i -b $annotation.promoter_4k_round_TSS.bed -wa -v | cut -f1-9,14,18-23 | awk '{OFS="\t"; print $0, "distal"}' >> $i.inpromoter
done

for i in *.commonuniq.inpromoter; do echo $i; cut -f4-5 $i | sort | uniq -c; done
#h3k4me3_peaks_intersected.Ct.bed.commonuniq.inpromoter
#   5010 common	distal
#  22001 common	proximal
#   2692 uniq	distal
#   3441 uniq	proximal
#h3k4me3_peaks_intersected.HD.bed.commonuniq.inpromoter
#   4411 common	distal
#  18956 common	proximal
#   1568 uniq	distal
#   2373 uniq	proximal
   
echo "Step 3: genes with H3k4kem3 peaks in proximal promoter"

# expression of genes with uniq/proximal peaks in each category

fgrep -f <(cat h3k4me3_peaks_intersected.HD.bed.commonuniq.inpromoter | grep uniq | grep proximal | intersectBed -a - -b $annotation.promoter_4k_round_TSS.bed -wao | sed 's/___/ /g' | cut -f2 -d' ' | uniq | sort -u) ../HD_DE_2014_03_04/no_filter/mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt > h3k4me3_peaks_intersected.HD.bed.uniq.proximal

fgrep -f <(cat h3k4me3_peaks_intersected.Ct.bed.commonuniq.inpromoter | grep uniq | grep proximal | intersectBed -a - -b $annotation.promoter_4k_round_TSS.bed -wao | sed 's/___/ /g' | cut -f2 -d' ' | uniq | sort -u) ../HD_DE_2014_03_04/no_filter/mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt > h3k4me3_peaks_intersected.Ct.bed.uniq.proximal &

R
df1=read.table("h3k4me3_peaks_intersected.Ct.bed.uniq.proximal", header=F)
df2=read.table("h3k4me3_peaks_intersected.HD.bed.uniq.proximal", header=F)
wilcox.test(df1$V3, df2$V3)
#
#	Wilcoxon rank sum test with continuity correction
#
#data:  df1$V3 and df2$V3
#W = 2990672, p-value = 3.108e-05
summary(df1$V3)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-2.18500 -0.31390 -0.06174 -0.02627  0.23000  3.52800 
summary(df2$V3)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-1.771000 -0.218200 -0.011970  0.007362  0.214700  2.382000 

## nearest genes with distal/uniq peaks in each category

bedtools slop -b 18000 -i $annotation.promoter_4k_round_TSS.bed -g hg19.genome > $annotation.promoter_2M_round_TSS.bed

# header:
# chr, start_peak, end_peak, uniq, distal, nearestGene, distance_TSS_to_mid_of_peak
cat h3k4me3_peaks_intersected.HD.bed.commonuniq.inpromoter | grep uniq | grep distal | intersectBed -a - -b $annotation.promoter_2M_round_TSS.bed -wo | awk '{d=$3+$2-$7-$8; if(d<0) d=-d; $12=int(d/2); OFS="\t"; print}' | sort -k1,1 -k2,2n -k3,3n -k12,12n | groupBy -g 1,2,3 -c 12 -o min -full |  cut -f1-5,9,12 > h3k4me3_peaks_intersected.HD.bed.uniq.distal.nearestGene.tab

cat h3k4me3_peaks_intersected.Ct.bed.commonuniq.inpromoter | grep uniq | grep distal | intersectBed -a - -b $annotation.promoter_2M_round_TSS.bed -wo | awk '{d=$3+$2-$7-$8; if(d<0) d=-d; $12=int(d/2); OFS="\t"; print}' | sort -k1,1 -k2,2n -k3,3n -k12,12n | groupBy -g 1,2,3 -c 12 -o min -full |  cut -f1-5,9,12 > h3k4me3_peaks_intersected.Ct.bed.uniq.distal.nearestGene.tab


############# (deprecated) ###########

# input
differentialPeaks=$1  #MAnorm_result.xls
peakfile1=$2  #$HOME/scratch/vs.input/h3k4me3_HD_M_peaks.xls
summit1=$3  #$HOME/scratch/vs.input/h3k4me3_HD_M_summits.bed
peakfile2=$4  #$HOME/scratch/vs.input/h3k4me3_control_M_g43yr_peaks.xls
summit2=$5  #$HOME/scratch/vs.input/h3k4me3_control_M_g43yr_summits.bed

significantPeaks=significantPeaks.bed

## filter peaks (FDR: 5% (then p-value is around 10^-50), fold_enrichment >4)
#grep -P "^chr[^\s]" $peakfile | paste - $summit |cut -f1-9,13,14 | awk '{if($9<10 && $8>4)print}' > $significantPeaks

echo "Step 1: get peaks from peak1 and peak2"
# use the differential peak calling
paste <(grep -P "^chr[^\s]" $peakfile1 | sort-bed -) <(sort-bed $summit1) <(fgrep peak1 $differentialPeaks | sort-bed -) > peak1.tab
paste <(grep -P "^chr[^\s]" $peakfile2 | sort-bed -) <(sort-bed $summit2) <(fgrep peak2 $differentialPeaks | sort-bed -) > peak2.tab
sort-bed peak1.tab peak2.tab > $significantPeaks

echo "Step 2: define peaks into proximal/distal"

echo -e "chr\tstart_peak\tend_peak\tlength\tsummit_pos\ttags_in_peak\tpvalue_MACS\tfold_enrichment\tFDR\tsummit_height\tpeak_type\traw_read_1\traw_read_2\tM_rescaled\tA_rescaled\tpvalue_MAnorm\tinpromoter" > $significantPeaks.inpromoter.tab
zcat $annotation | fgrep -w transcript | sed 's/[";]//g;' | awk '{OFS="\t"; tss=($7=="+")?($4-1):($5-1); if(tss<2000) tss=2000; print $1, tss-2000,tss+2000,$12"___"$10"___"$14"___"$18}' | bedtools intersect -b stdin -a $significantPeaks -wa -u | cut -f1-9,14,18-23 | awk '{OFS="\t"; print $0, "proximal"}' >> $significantPeaks.inpromoter.tab
zcat $annotation | fgrep -w transcript | sed 's/[";]//g;' | awk '{OFS="\t"; tss=($7=="+")?($4-1):($5-1); if(tss<2000) tss=2000; print $1, tss-2000,tss+2000,$12"___"$10"___"$14"___"$18}' | bedtools intersect -b stdin -a $significantPeaks -wa -v | cut -f1-9,14,18-23 | awk '{OFS="\t"; print $0, "distal"}' >> $significantPeaks.inpromoter.tab

echo "Step 3: genes with H3k4kem3 peaks in proximal promoter"

## intersect with [-2k, 2k] region of annotated TSS/gene (Caution: the TSS in the gene model in Gencode is actually the most-upstream TSS among all isoforms)
#zcat $annotation | fgrep -w gene | sed 's/[";]//g;' | awk '{OFS="\t"; tss=($7=="+")?($4-1):($5-1); if(tss<2000) tss=2000; print $1, tss-2000,tss+2000,$10"___"$14"___"$18}' | bedtools intersect -a stdin -b <(awk '{OFS="\t"; $13=$18; print}' $significantPeaks | cut -f1-9,13,14) -wao | awk '{OFS="\t"; dis_tss2sumit=$6+$9-($2+$3)/2; if(dis_tss2sumit<0) dis_tss2sumit=-dis_tss2sumit; print $0, dis_tss2sumit}' | sort -k4,4 -k17,17n | awk '{if($4!=id) {id=$4; print;}}' | sed 's/___/\t/g' > $HOME/projects/bu_neuro/data/genes.H3k4me3inpromoter.tab
#zcat $annotation | fgrep -w transcript | sed 's/[";]//g;' | awk '{OFS="\t"; tss=($7=="+")?($4-1):($5-1); if(tss<2000) tss=2000; print $1, tss-2000,tss+2000,$12"___"$10"___"$14"___"$18}' | bedtools intersect -a stdin -b <(awk '{OFS="\t"; $13=$18; print}' $significantPeaks | cut -f1-9,13,14) -wao | awk '{OFS="\t"; dis_tss2sumit=$6+$9-($2+$3)/2; if(dis_tss2sumit<0) dis_tss2sumit=-dis_tss2sumit; print $0, dis_tss2sumit}' | sort -k4,4 -k17,17n | awk '{if($4!=id) {id=$4; print;}}' | sed 's/___/\t/g' > $HOME/projects/bu_neuro/data/transcripts.H3k4me3inpromoter.tab

## assign the closest peak to each TSS first, and then assigning the overall gene expression to the strongest peak in ALL isoforms' promoter.
zcat $annotation | fgrep -w transcript | sed 's/[";]//g;' | awk '{OFS="\t"; tss=($7=="+")?($4-1):($5-1); if(tss<2000) tss=2000; print $1, tss-2000,tss+2000,$12"___"$10"___"$14"___"$18}' | bedtools intersect -a stdin -b $significantPeaks -wao | cut -f1-13,18,22-27 | awk '{OFS="\t"; dis_tss2sumit=$6+$9-($2+$3)/2; if(dis_tss2sumit<0) dis_tss2sumit=-dis_tss2sumit; print $0, dis_tss2sumit}' | sort -k4,4 -k21,21n | awk '{if($4!=id) {id=$4; print;}}' | sed 's/___/\t/g' | sort -k5,5 -k17,17nr | awk '{if($5!=id) {id=$5; print;}}' > genes.H3k4me3inpromoter.tab
#zcat $annotation | fgrep -w transcript | sed 's/[";]//g;' | awk '{OFS="\t"; tss=($7=="+")?($4-1):($5-1); if(tss<2000) tss=2000; print $1, tss-2000,tss+2000,$12"___"$10"___"$14"___"$18}' | bedtools intersect -a stdin -b $significantPeaks -wao | awk '{OFS="\t"; dis_tss2sumit=$6+$9-($2+$3)/2; if(dis_tss2sumit<0) dis_tss2sumit=-dis_tss2sumit; print $0, dis_tss2sumit}' | sort -k4,4 -k29,29n | awk '{if($4!=id) {id=$4; print;}}' | sed 's/___/\t/g' > $HOME/projects/bu_neuro/data/transcripts.H3k4me3inpromoter.tab

echo "Step 4: merge with expression level"
# merge with expression level

#echo -e "ID\tchr\tstart_promoter\tend_promoter\ttype\tname\tstart_peak\tend_peak\tlength\tsummit_pos\ttags_in_peak","-10log10(pvalue)","fold_enrichment\tFDR(%)\tpeak_type\tsummit_height\tC_0014\tC_0018\tC_0021\tC_0029\tC_0031\tC_0032\tC_0033\tC_0034\tC_0035\tC_0036\tC_0037\tC_0038\tC_0039\tC_0050\tC_0053\tC_0070\tC_0077\tC_0087\tH_0001\tH_0002\tH_0003\tH_0005\tH_0006\tH_0007\tH_0008\tH_0009\tH_0010\tH_0012\tH_0013\tH_0014\tH_0539\tH_0657\tH_0658\tH_0681\tH_0695\tH_0700\tH_0726\tH_0740\tH_0750" > $HOME/projects/bu_neuro/data/genes.H3k4me3inpromoter.expression.tab
join -1 5 -2 1 -a 1 genes.H3k4me3inpromoter.tab <(sed 's/"//g' $expression | sort -k1,1) | sed 's/ /\t/g' > allgenes.H3k4me3.expression.v2
