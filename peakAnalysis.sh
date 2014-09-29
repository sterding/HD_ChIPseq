#!/bin/bash

cd /Users/xdong/projects/HD/data/scratch
# required files
annotation=../gencode.v17.annotation.gtf.gz
#gzcat $annotation | fgrep -w transcript | sed 's/[";]//g;' | awk '{OFS="\t"; tss=($7=="+")?($4-1):($5-1); if(tss<2000) tss=2000; print $1, tss-2000,tss+2000,$12"___"$10"___"$14"___"$18, 4000, $7}'  > $annotation.promoter_4k_round_TSS.bed
# bedtools slop -b 998000 -i $annotation.promoter_4k_round_TSS.bed -g hg19.chrom.sizes > $annotation.promoter_1M_round_TSS.bed

#expression=$HOME/projects/bu_neuro/data/mlhd_DESeq_normalized_counts_all_controls.txt  # V1 (21HD and 54Ct)
#expression=$HOME/projects/bu_neuro/data/mlhd_DESeq_normalized_counts_all_controls_strict.txt  # V2 (21HD and 50Ct, and using strict DEseq criteria)

# input
peakfile1=h3k4me3_peaks_intersected.Ct.nofilter.bed 
peakfile2=h3k4me3_peaks_intersected.HD.nofilter.bed

echo "Step 0: merge peaks into a union set and then quantify their signal"
cat  $peakfile1 $peakfile2 | awk '$4>=3' |  sortBed | mergeBed | awk '($3-$2)>100' > h3k4me3_peaks_union.bed

intersectBed -a h3k4me3_peaks_union.bed -b h3k4me3_peaks_intersected.Ct.nofilter.bed -wao | sed 's/\./0/g' | groupBy -g 1,2,3 -c 7 -o max > h3k4me3_peaks_union.Ct_depth.txt
intersectBed -a h3k4me3_peaks_union.bed -b h3k4me3_peaks_intersected.HD.nofilter.bed -wao | sed 's/\./0/g' | groupBy -g 1,2,3 -c 7 -o max > h3k4me3_peaks_union.HD_depth.txt

## there is a BUG in bigwigsummary!!
#../../src/toBinRegionsOnBigwig.sh trimmedmean.normalized.Ct.bw h3k4me3_peaks_union.bed 1 max > h3k4me3_peaks_union.Ct_signal.txt
#../../src/toBinRegionsOnBigwig.sh trimmedmean.normalized.HD.bw h3k4me3_peaks_union.bed 1 max > h3k4me3_peaks_union.HD_signal.txt
#paste h3k4me3_peaks_union.bed h3k4me3_peaks_union.Ct_signal.txt h3k4me3_peaks_union.HD_signal.txt | sed 's/ /\t/g' | cut -f1-3,5,7 | sortBed > h3k4me3_peaks_union.signal.bed

# change to use intersectBed + groupBy
intersectBed -a h3k4me3_peaks_union.bed -b trimmedmean.normalized.Ct.bedGraph -wao -sorted | sed 's/\t\.\t/\t0\t/g' | groupBy -g 1,2,3 -c 7 -o max > h3k4me3_peaks_union.Ct_signal.txt
intersectBed -a h3k4me3_peaks_union.bed -b trimmedmean.normalized.HD.bedGraph -wao -sorted | sed 's/\t\.\t/\t0\t/g' | groupBy -g 1,2,3 -c 7 -o max > h3k4me3_peaks_union.HD_signal.txt
paste h3k4me3_peaks_union.Ct_signal.txt h3k4me3_peaks_union.HD_signal.txt | sed 's/ /\t/g' | cut -f1-4,8 | sortBed | awk '{OFS="\t"; $4=$4+0; $5=$5+0;print}'> h3k4me3_peaks_union.signal.bed

# update: if peakCenter is upstream of TSS, then the distance is minus; otherwise, plus
bedtools intersect -a h3k4me3_peaks_union.bed -b $annotation.promoter_1M_round_TSS.bed -wao | awk '{OFS="\t"; tss=(($6-$5)==2e6)?(($5+$6)/2):(($5==0)?($6-1e6):($5+1e6)); d=tss-($3+$2)/2; if($9=="-") d=-d; d2=d;if(d2<0)d2=-d2;print $0,d,d2;}' | sort -k1,1 -k2,2n -k12,12n | awk '{split($7, a, "___");b=a[2]"___"a[3]"___"a[4]; OFS="\t";if(id!=$1$2$3) {if($4!=".") print $1, $2, $3,b,$11; else print $1,$2,$3,".","."; id=$1$2$3;} }' > h3k4me3_peaks_union.neighborhood.bed

paste h3k4me3_peaks_union.signal.bed <(cut -f4 h3k4me3_peaks_union.Ct_depth.txt) <(cut -f4 h3k4me3_peaks_union.HD_depth.txt) <(cut -f4-5 h3k4me3_peaks_union.neighborhood.bed) > h3k4me3_peaks_union.signal.neighborhood.bed

echo -e "chr\tstart_peak\tend_peak\trpm_summit_CT\trpm_summit_HD\tsamples_support_peak_CT\tsamples_support_peak_HD\tnearest_gene\tdistance_peakcenter2TSS\trelative_location" > h3k4me3_peaks_union.signal.neighborhood.annotated.bed.xls

gzcat $annotation | fgrep -w gene | sed 's/[";]//g;' | cut -f1,4,5 | intersectBed -a h3k4me3_peaks_union.signal.neighborhood.bed -b - -wao | cut -f1-10| uniq | awk '{OFS="\t"; $10=($10==".")?"intergenic":"intragenic"; print;}' >> h3k4me3_peaks_union.signal.neighborhood.annotated.bed.xls

cp h3k4me3_peaks_union.signal.neighborhood.annotated.bed.xls ~/Dropbox/huntington_bu/data/


#############################################################
# how much RNA-defined eRNA overlap with other enhancers --> venn diagram
#############################################################
grep -w distal ~/Dropbox/huntington_bu/data/h3k4me3_peaks_union.signal.neighborhood.annotated.bed.xls | awk '{OFS="\t"; print $2,$3,$4,$1}' > distalPeaks.bed

# Roadmap Epigenomics enhancers for Brain Dorsolateral Prefrontal Cortex (https://sites.google.com/site/anshulkundaje/projects/epigenomeroadmap)
curl -s http://www.broadinstitute.org/~anshul/projects/roadmap/segmentations/models/coreMarks/parallel/set2/final/E073_15_coreMarks_segments.bed | awk '$4~/E6|E7|E12/' | intersectBed -a distalPeaks.bed -b stdin -c | sort -k4,4 > distalPeaks.overlap.txt

## super enhancer in brain
intersectBed -a distalPeaks.bed -b ../enhancer.superenhancer.hg19.BI_Brain_Mid_Frontal_Lobe.bed -c | sort -k4,4 | cut -f5 | paste distalPeaks.overlap.txt - > tmp.list
mv tmp.list distalPeaks.overlap.txt

# CAGE-defined enhancers
curl -s http://enhancer.binf.ku.dk/presets/permissive_enhancers.bed | fgrep -v track | intersectBed -a distalPeaks.bed -b stdin -c | sort -k4,4 | cut -f5 | paste distalPeaks.overlap.txt - > tmp.list
mv tmp.list distalPeaks.overlap.txt

# CAGE-defined enhancer expressed in 
enhancers_CAGE=~/projects/HD/data/hg19_permissive_enhancers_expression_rle_tpm.brain.tab  # in total, 43011 enhancers defined by CAGE
# wget http://enhancer.binf.ku.dk/presets/hg19_permissive_enhancers_expression_rle_tpm.csv.gz
# rowsToCols hg19_permissive_enhancers_expression_rle_tpm.csv -fs=, stdout | grep -E "chr1:858256-858648|CNhs10617|CNhs10647" | rowsToCols -tab stdin stdout | sed 's/"//g' | awk '{OFS="\t"; print $3, $1,$2;}' | sed 's/[:-]/ \t/g' | grep chr > hg19_permissive_enhancers_expression_rle_tpm.brain.tab
# add header "#chr	start	end	TPM_frontallobeadult_CNhs10647	TPM_brainadult_CNhs10617"
awk '{OFS="\t"; if($4>0) print $1,$2,$3,$1"_"$2"_"$3"_"$4;}' $enhancers_CAGE | intersectBed -a distalPeaks.bed -b stdin -c | sort -k4,4 | cut -f5 | paste distalPeaks.overlap.txt - > tmp.list
mv tmp.list distalPeaks.overlap.txt

## DNase cluster
#intersectBed -a distalPeaks.bed -b ../DNase/DNase.distal.bed -c | sort -k4,4 | cut -f5 | paste distalPeaks.overlap.txt - > tmp.list
#mv tmp.list distalPeaks.overlap.txt
#
## TFBS
#intersectBed -a distalPeaks.bed -b ../TFBS/TFBS.distal.bed -c | sort -k4,4 | cut -f5 | paste distalPeaks.overlap.txt - > tmp.list
#mv tmp.list distalPeaks.overlap.txt
#
## Conservation
#intersectBed -a distalPeaks.bed -b ../Conservation/Conservation.distal.bed -c | sort -k4,4 | cut -f5 | paste distalPeaks.overlap.txt - > tmp.list
#mv tmp.list distalPeaks.overlap.txt

## overlap %
cat distalPeaks.overlap.txt | datamash mean 5 mean 6 mean 7 mean 8 
#0.21051382461137	0.15343538625845	0.17408330366679	0.12175151299395

# fisher exact test to see the signifiance
# generate a random background region with the same length distribution
../../src/toGenerateRandomRegions.sh distalPeaks.bed | awk '{OFS="\t"; print $1, $2,$3,$1"_"$2"_"$3}'> randomPeaks.bed

curl -s http://www.broadinstitute.org/~anshul/projects/roadmap/segmentations/models/coreMarks/parallel/set2/final/E073_15_coreMarks_segments.bed | awk '$4~/E6|E7|E12/' | intersectBed -a randomPeaks.bed -b stdin -c | sort -k4,4 > randomPeaks.overlap.txt
intersectBed -a randomPeaks.bed -b ../enhancer.superenhancer.hg19.BI_Brain_Mid_Frontal_Lobe.bed -c | sort -k4,4 | cut -f5 | paste randomPeaks.overlap.txt - > tmp.list
mv tmp.list randomPeaks.overlap.txt
curl -s http://enhancer.binf.ku.dk/presets/permissive_enhancers.bed | fgrep -v track | intersectBed -a randomPeaks.bed -b stdin -c | sort -k4,4 | cut -f5 | paste randomPeaks.overlap.txt - > tmp.list
mv tmp.list randomPeaks.overlap.txt
awk '{OFS="\t"; if($4>0) print $1,$2,$3,$1"_"$2"_"$3"_"$4;}' $enhancers_CAGE | intersectBed -a randomPeaks.bed -b stdin -c | sort -k4,4 | cut -f5 | paste randomPeaks.overlap.txt - > tmp.list
mv tmp.list randomPeaks.overlap.txt
cat randomPeaks.overlap.txt | datamash mean 5 mean 6 mean 7 mean 8 
# 0.045923816304735	0.0084252996321348	0.015189272576243	0.0034413195680551




######## DEPRECATED #############

echo "Step 1: classify peaks into common/uniq"

intersectBed -a $peakfile1 -b $peakfile2 -wao | awk '{OFS="\t"; print $0, ($5==".")?"uniq":"common";}' > $peakfile1.commonuniq
intersectBed -a $peakfile2 -b $peakfile1 -wao | awk '{OFS="\t"; print $0, ($5==".")?"uniq":"common";}' > $peakfile2.commonuniq

for i in *commonuniq; do echo $i; cut -f1-3,10 $i | sort -u | cut -f4 | sort | uniq -c; done
#h3k4me3_peaks_intersected.Ct.bed.commonuniq
#  28182 common
#   4962 uniq
#h3k4me3_peaks_intersected.HD.bed.commonuniq
#  24775 common
#   2533 uniq
   
echo "Step 2: classify peaks into proximal/distal"

for i in *commonuniq; do
    echo $i; 
    bedtools intersect -a $i -b $annotation.promoter_4k_round_TSS.bed -wa -u | cut -f1-9,14,18-23 | awk '{OFS="\t"; print $0":proximal"}' > $i.inpromoter
    bedtools intersect -a $i -b $annotation.promoter_4k_round_TSS.bed -wa -v | cut -f1-9,14,18-23 | awk '{OFS="\t"; print $0":distal"}' >> $i.inpromoter
done

for i in *.commonuniq.inpromoter; do echo $i; cut -f4-5 $i | sort | uniq -c; done
#h3k4me3_peaks_intersected.Ct.bed.commonuniq.inpromoter
#   5798 common	distal
#  24368 common	proximal
#   2286 uniq	distal
#   2676 uniq	proximal
#h3k4me3_peaks_intersected.HD.bed.commonuniq.inpromoter
#   5517 common	distal
#  24649 common	proximal
#   1136 uniq	distal
#   1397 uniq	proximal
   
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
#W = 1316863, p-value = 4.303e-06
summary(df1$V3)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-2.18500 -0.31570 -0.06465 -0.02600  0.23120  3.52800 
summary(df2$V3)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-1.330000 -0.204300  0.001233  0.023860  0.238000  2.382000 

## nearest genes with distal/uniq peaks in each category
## --------------------------------------------

bedtools slop -b 18000 -i $annotation.promoter_4k_round_TSS.bed -g hg19.genome > $annotation.promoter_2M_round_TSS.bed

# header:
# chr, start_peak, end_peak, uniq, distal, nearestGene, distance_TSS_to_mid_of_peak
cat h3k4me3_peaks_intersected.HD.bed.commonuniq.inpromoter | grep uniq | grep distal | intersectBed -a - -b $annotation.promoter_2M_round_TSS.bed -wo | awk '{d=$3+$2-$7-$8; if(d<0) d=-d; $12=int(d/2); OFS="\t"; print}' | sort -k1,1 -k2,2n -k3,3n -k12,12n | groupBy -g 1,2,3 -c 12 -o min -full |  cut -f1-5,9,12 > h3k4me3_peaks_intersected.HD.bed.uniq.distal.nearestGene.tab

cat h3k4me3_peaks_intersected.Ct.bed.commonuniq.inpromoter | grep uniq | grep distal | intersectBed -a - -b $annotation.promoter_2M_round_TSS.bed -wo | awk '{d=$3+$2-$7-$8; if(d<0) d=-d; $12=int(d/2); OFS="\t"; print}' | sort -k1,1 -k2,2n -k3,3n -k12,12n | groupBy -g 1,2,3 -c 12 -o min -full |  cut -f1-5,9,12 > h3k4me3_peaks_intersected.Ct.bed.uniq.distal.nearestGene.tab

## nearest genes with all distal peaks (TO CHANGE)
## --------------------------------------------
cat *.commonuniq.inpromoter | grep distal | intersectBed -a - -b $annotation.promoter_2M_round_TSS.bed -wo | awk '{d=$3+$2-$7-$8; if(d<0) d=-d; $12=int(d/2); OFS="\t"; print}' | sort -k1,1 -k2,2n -k3,3n -k12,12n | groupBy -g 1,2,3 -c 12 -o min -full |  cut -f1-5,9,12 > h3k4me3_peaks_intersected.distal.nearestGene.tab


## TF (from ENCODE) occurnaces per distal peaks
## --------------------------------------------
# clustered TFBS (count all Peaks for 161 transcription factors in 91 cell types)
#wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz
#gunzip -c wgEncodeRegTfbsClusteredV3.bed.gz | awk '{OFS="\t"; print $1,$2,$3,$4,0,"+",$2,$3,"0",$6,$7,$8}' > wgEncodeRegTfbsClusteredV3.bed12

bed12ToBed6 -i ../wgEncodeRegTfbsClusteredV3.bed12 | bedtools intersect -a - -b $annotation.promoter_4k_round_TSS.bed -v > ../wgEncodeRegTfbsClusteredV3.bed12.distal
grep "common:distal" h3k4me3_peaks_intersected.HD.bed.commonuniq.inpromoter | intersectBed -a ../wgEncodeRegTfbsClusteredV3.bed12.distal -b - -u 

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
