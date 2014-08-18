#################### stratedgy 1: merge peaks called from individual alignment ###########
cd ~/scratch

# Note: the 6 h3k4me3 samples in batch1 have better quality
# a peak is defined as a region which was called as peak at least 50% samples (e.g. 3 out of 6 for this case)
multiIntersectBed -i h3k4me3_HD_batch1_*/*peaks.bed | awk '{OFS="\t";if($4>=3)print $1,$2,$3,".",$4}' | mergeBed -scores max > h3k4me3_peaks_intersected.HD.bed

# same for control (3 out 6)
multiIntersectBed -i Schraham_1644/*peaks.bed Schraham_1713/*peaks.bed Schraham_3706R/*peaks.bed Schraham_7244/*peaks.bed Schraham_R30/*peaks.bed Schraham_R7/*peaks.bed | awk '{OFS="\t";if($4>=3)print $1,$2,$3,".",$4}' | mergeBed -scores max > h3k4me3_peaks_intersected.Ct.bed

## also record all union peaks (e.g. peak called in at least one sample)
multiIntersectBed -i h3k4me3_HD_batch1_*/*peaks.bed | cut -f1-5 > h3k4me3_peaks_intersected.HD.nofilter.bed
multiIntersectBed -i Schraham_*/*peaks.bed | cut -f1-5  > h3k4me3_peaks_intersected.Ct.nofilter.bed


##################### stratedgy 2: peak calling for merged samples (vs. 2037input) #####
## merge alignment
#cat Schraham_{R7,3706R,7244,R30,1644}/accepted_hits.bed > h3k4me3_control_M_g43yr.accepted_hits.bed
#ln -s h3k4me3_control_M_g43yr.accepted_hits.bed h3k4me3_Ct_Schraham.accepted_hits.bed # only 5 Male and >43yr samples
#cat h3k4me3_HD_batch1_*/accepted_hits.bed > h3k4me3_HD_batch1.accepted_hits.bed
#
#controlIP=Schraham_2037in/accepted_hits.bed
#macs14 -t h3k4me3_Ct_Schraham.accepted_hits.bed -c $controlIP -f BED -n h3k4me3_Ct_Schraham 2> h3k4me3_Ct_Schraham.macs.log &
#macs14 -t h3k4me3_HD_batch1.accepted_hits.bed -c $controlIP -f BED -n h3k4me3_HD_batch1 2> h3k4me3_HD_batch1.macs.log &
#
######  calling differential binding between samples using MAnorm #####
## HD vs. Ct  eaks_intersected.bed ../h3k4me3_Ct_Schraham_peaks_intersected.bed ../h3k4me3_HD_batch1.accepted_hits.bed ../h3k4me3_control_M_g43yr.accepted_hits.bed 140 140 &
#
#$HOME/projects/bu_neuro/src/peakAnalysis2.sh MAnorm_result.xls
#mkdir pdf
#Rscript $HOME/projects/bu_neuro/src/peakAnalysis2.R $MAnorm_result.inpromoter.tab
#Rscript $HOME/projects/bu_neuro/src/geneAnalysis2.R
#mv *.pdf pdf

#################### normalized bigwig and merge samples ###########
for i in Schraham_* h3k4me3_HD_batch1_*; do
    n=`wc -l $i/accepted_hits.bed| cut -f1 -d' '`
    awk -vn=$n 'BEGIN{OFS="\t"; print "#total_reads="n;}{$4=$4*1e6/n; print}' $i/accepted_hits.bedGraph > $i/accepted_hits.normalized.bedGraph
    bedGraphToBigWig $i/accepted_hits.normalized.bedGraph hg19.chrom.sizes $i/accepted_hits.normalized.bw
    ln -sf $i/accepted_hits.normalized.bw $i.accepted_hits.normalized.bw
done

unionBedGraphs -i `ls Schraham_*/*normalized.bedGraph` | awk 'function trimmedMean(v, p) {c=asort(v,j); a=int(c*p); s=0; for(i=a+1;i<=(c-a);i++) s+=j[i];return s/(c-2*a);} {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; tm=trimmedMean(S, 0.1); if(tm!=0) print $1,$2,$3,tm}' > trimmedmean.normalized.Ct.bedGraph && bedGraphToBigWig trimmedmean.normalized.Ct.bedGraph hg19.chrom.sizes trimmedmean.normalized.Ct.bw

unionBedGraphs -i `ls h3k4me3_HD_batch1_*/*normalized.bedGraph` | awk 'function trimmedMean(v, p) {c=asort(v,j); a=int(c*p); s=0; for(i=a+1;i<=(c-a);i++) s+=j[i];return s/(c-2*a);} {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; tm=trimmedMean(S, 0.1); if(tm!=0) print $1,$2,$3,tm}' > trimmedmean.normalized.HD.bedGraph && bedGraphToBigWig trimmedmean.normalized.HD.bedGraph hg19.chrom.sizes trimmedmean.normalized.HD.bw

scp *normalized.bw dongx@zlab.umassmed.edu:/home/dongx/public_html/tracks/HD/h3k4me3/merged


#################### peak calling for merged alignment (deprecated) ###########

cd ~/scratch

##### merge controls #####

# merge control (M, >43yr)
cat Schraham_{R7,3706R,7244,R30,1644}/accepted_hits.bed > h3k4me3_control_M_g43yr.accepted_hits.bed
# merge control (M, <43yr)
cat Schraham_{j13,i14,k6,1698}/accepted_hits.bed > h3k4me3_control_M_l43yr.accepted_hits.bed
# merge control (M)
cat h3k4me3_control_M_g43yr.accepted_hits.bed h3k4me3_control_M_l43yr.accepted_hits.bed > h3k4me3_control_M.accepted_hits.bed
# merge control (F, >43yr)
cat Schraham_{1713,1796,2037,3520R}/accepted_hits.bed > h3k4me3_control_F_g43yr.accepted_hits.bed
# merge control (F, <43yr)
ln -s Schraham_B24/accepted_hits.bed h3k4me3_control_F_l43yr.accepted_hits.bed
# merge control (F)
cat h3k4me3_control_F_g43yr.accepted_hits.bed Schraham_B24/accepted_hits.bed > h3k4me3_control_F.accepted_hits.bed
ln -s Schraham_1796/accepted_hits.bed h3k4me3_control_neuro-.accepted_hits.bed
cat Schraham_*/accepted_hits.bed > h3k4me3_control_combined.accepted_hits.bed

# merge H3k4me3 neuro+
cat `ls combined_Lane*/accepted_hits.bed h3k4me3_HD_batch1_*/*hits.bed | grep -v 3584` > h3k4me3_HD_M.accepted_hits.bed
ln -s h3k4me3_HD_batch1_3584/accepted_hits.bed h3k4me3_HD_F.accepted_hits.bed

##### peak calling (vs. 2037input) #####
controlIP=Schraham_2037in/accepted_hits.bed
macs14 -t h3k4me3_HD_M.accepted_hits.bed -c $controlIP -f BED -n h3k4me3_HD_M 2> h3k4me3_HD_M_macs.log &
macs14 -t h3k4me3_HD_F.accepted_hits.bed -c $controlIP -f BED -n h3k4me3_HD_F 2> h3k4me3_HD_F_macs.log &
macs14 -t h3k4me3_control_M_g43yr.accepted_hits.bed -c $controlIP -f BED -n h3k4me3_control_M_g43yr 2> h3k4me3_control_M_g43yr_macs.log &
macs14 -t h3k4me3_control_M_l43yr.accepted_hits.bed -c $controlIP -f BED -n h3k4me3_control_M_l43yr 2> h3k4me3_control_M_l43yr_macs.log &
macs14 -t h3k4me3_control_M.accepted_hits.bed -c $controlIP -f BED -n h3k4me3_control_M 2> h3k4me3_control_M_macs.log &
macs14 -t h3k4me3_control_F_g43yr.accepted_hits.bed -c $controlIP -f BED -n h3k4me3_control_F_g43yr 2> h3k4me3_control_F_g43yr_macs.log &
macs14 -t h3k4me3_control_F_l43yr.accepted_hits.bed -c $controlIP -f BED -n h3k4me3_control_F_l43yr 2> h3k4me3_control_F_l43yr_macs.log &
macs14 -t h3k4me3_control_F.accepted_hits.bed -c $controlIP -f BED -n h3k4me3_control_F 2> h3k4me3_control_F_macs.log &
mkdir vs.input;
find h3k4me3_HD_M_* h3k4me3_HD_F_* h3k4me3_control_M_g43yr_* h3k4me3_control_M_l43yr_* h3k4me3_control_M_* h3k4me3_control_F_g43yr_* h3k4me3_control_F_l43yr_* h3k4me3_control_F_* -not -name \*accepted_hits.bed -exec mv {} vs.input/ \;

#####  differential binding between samples #####
# HD vs. Ct
mkdir ~/scratch/HD.vs.Ct_M_g43yr
cd ~/scratch/HD.vs.Ct_M_g43yr
export MAnorm=$HOME/bin/MAnorm_Linux_R_Package/
$MAnorm/MAnorm2.sh ../vs.input/h3k4me3_HD_M_peaks.bed ../vs.input/h3k4me3_control_M_g43yr_peaks.bed ../h3k4me3_HD_M.accepted_hits.bed ../h3k4me3_control_M_g43yr.accepted_hits.bed 140 108 &
$HOME/projects/bu_neuro/src/peakAnalysis.sh MAnorm_result.xls ../vs.input/h3k4me3_HD_M_peaks.xls ../vs.input/h3k4me3_HD_M_summits.bed ../vs.input/h3k4me3_control_M_g43yr_peaks.xls ../vs.input/h3k4me3_control_M_g43yr_summits.bed
Rscript $HOME/projects/bu_neuro/src/peakAnalysis.R
Rscript $HOME/projects/bu_neuro/src/geneAnalysis.R


mkdir ~/scratch/HD.vs.Ct
cd ~/scratch/HD.vs.Ct
export MAnorm=$HOME/bin/MAnorm_Linux_R_Package/
$MAnorm/MAnorm2.sh ../vs.input/h3k4me3_HD_M_peaks.bed ../vs.input/h3k4me3_control_M_peaks.bed ../h3k4me3_HD_M.accepted_hits.bed ../h3k4me3_control_M.accepted_hits.bed 140 106 &
$HOME/projects/bu_neuro/src/peakAnalysis.sh MAnorm_result.xls ../vs.input/h3k4me3_HD_M_peaks.xls ../vs.input/h3k4me3_HD_M_summits.bed ../vs.input/h3k4me3_control_M_peaks.xls ../vs.input/h3k4me3_control_M_summits.bed
Rscript $HOME/projects/bu_neuro/src/peakAnalysis.R
Rscript $HOME/projects/bu_neuro/src/geneAnalysis.R

mkdir  ~/scratch/HD.vs.Ct_M_l43yr
cd ~/scratch/HD.vs.Ct_M_l43yr
export MAnorm=$HOME/bin/MAnorm_Linux_R_Package/
$MAnorm/MAnorm2.sh ../vs.input/h3k4me3_HD_M_peaks.bed ../vs.input/h3k4me3_control_M_l43yr_peaks.bed ../h3k4me3_HD_M.accepted_hits.bed ../h3k4me3_control_M_l43yr.accepted_hits.bed 140 200 &


#################### tranasfer data/track to zlab ###########

cd ~/scratch
for i in combined_Lane*/*bam* combined_Lane*/*bw combined_Lane*/*summary; do scp $i zlab:~/public_html/tracks/HD/h3k4me3/Rick/${i/\//_}; done
for i in combined_Lane*/*summits.bed combined_Lane*/*peaks.bed combined_Lane*/*xls; do scp $i zlab:~/public_html/tracks/HD/h3k4me3/Rick/; done
for i in h3k4me3_HD_batch1_*/*bam* h3k4me3_HD_batch1_*/*bw h3k4me3_HD_batch1_*/*summary; do scp $i zlab:~/public_html/tracks/HD/h3k4me3/JF/${i/\//_}; done
for i in h3k4me3_HD_batch1_*/*xls h3k4me3_HD_batch1_*/*summits.bed h3k4me3_HD_batch1_*/*peaks.bed; do scp $i zlab:~/public_html/tracks/HD/h3k4me3/JF/; done
for i in JF_h3k4me3_*/*bam* JF_h3k4me3_*/*bw JF_h3k4me3_*/*summary; do scp $i zlab:~/public_html/tracks/HD/h3k4me3/JF/${i/\//_}; done
for i in JF_h3k4me3_*/*xls JF_h3k4me3_*/*summits.bed JF_h3k4me3_*/*peaks.bed; do scp $i zlab:~/public_html/tracks/HD/h3k4me3/JF/; done
for i in Schraham_*/*bam* Schraham_*/*bw Schraham_*/*summary; do scp $i zlab:~/public_html/tracks/HD/h3k4me3/Schraham/${i/\//_}; done
for i in Schraham_*/*xls Schraham_*/*summits.bed Schraham_*/*peaks.bed; do scp $i zlab:~/public_html/tracks/HD/h3k4me3/Schraham/; done


> $HOME/projects/bu_neuro/data/trackDb.ChIPseq.txt
for i in ~/nearline/BU/h3k4me3/chipseq_h3k4me3_May2013/130420-376_and_130429-380_Richard_Myers/*
do
    j=`basename $i`;
    samplename=${j/.*/}
    index="hg19"; color="0,0,255"
    ### generate trackDb.txt
    echo "
    ######################################################
    ##  $samplename
    ######################################################
    ## ---------- BAM
    #track $samplename.bam
    #bigDataUrl http://zlab.umassmed.edu/~dongx/tracks/HD/h3k4me3/Rick/${samplename}_accepted_hits.bam
    #shortLabel $samplename.bam
    #longLabel $samplename.bam
    #type bam
    #pairEndsByName .
    #pairSearchRange 10000
    #bamColorMode strand
    #maxWindowToDraw 200000
    #db $index
    #visibility pack
    #colorByStrand 200,100,0 0,100,200

    ## ---------- bigwig
    track $samplename.bw
    bigDataUrl http://zlab.umassmed.edu/~dongx/tracks/HD/h3k4me3/Rick/${samplename}_accepted_hits.bw
    shortLabel $samplename.bw
    longLabel $samplename.bw
    type bigWig
    autoScale on
    maxHeightPixels 100:20:8
    visibility full
    aggregate transparentOverlay
    showSubtrackColorOnUi on
    windowingFunction mean+whiskers
    priority 1.4
    configurable on
    color $color

    " | sed 's/^[ \t]*//g' >> $HOME/projects/bu_neuro/data/trackDb.ChIPseq.txt
done

for i in ~/nearline/BU/h3k4me3/chipseq_JF/*.gz ~/nearline/BU/h3k4me3/h3k4me3_HD_batch1/*.fq
do
    j=`basename $i`;
    samplename=${j/.*/}
    [[ $i =~ /h3k4me3_HD_batch1/ ]] && samplename='h3k4me3_HD_batch1_'$samplename
    index="hg19"; color="0,0,100"
    ### generate trackDb.txt
    echo "
    ######################################################
    ##  $samplename
    ######################################################
    ## ---------- BAM
    #track $samplename.bam
    #bigDataUrl http://zlab.umassmed.edu/~dongx/tracks/HD/h3k4me3/JF/${samplename}_accepted_hits.bam
    #shortLabel $samplename.bam
    #longLabel $samplename.bam
    #type bam
    #pairEndsByName .
    #pairSearchRange 10000
    #bamColorMode strand
    #maxWindowToDraw 200000
    #db $index
    #visibility pack
    #colorByStrand 200,100,0 0,100,200

    ## ---------- bigwig
    track $samplename.bw
    bigDataUrl http://zlab.umassmed.edu/~dongx/tracks/HD/h3k4me3/JF/${samplename}_accepted_hits.bw
    shortLabel $samplename.bw
    longLabel $samplename.bw
    type bigWig
    autoScale on
    maxHeightPixels 100:20:8
    visibility full
    aggregate transparentOverlay
    showSubtrackColorOnUi on
    windowingFunction mean+whiskers
    priority 1.4
    configurable on
    color $color

    " | sed 's/^[ \t]*//g' >> $HOME/projects/bu_neuro/data/trackDb.ChIPseq.txt
done

for i in ~/nearline/BU/h3k4me3/control_from_Schraham/*.fq
do
    j=`basename $i`;
    samplename=${j/.*/}
    samplename='Schraham_'$samplename
    index="hg19"; color="100,100,100"
    ### generate trackDb.txt
    echo "
    ######################################################
    ##  $samplename
    ######################################################
    ## ---------- BAM
    #track $samplename.bam
    #bigDataUrl http://zlab.umassmed.edu/~dongx/tracks/HD/h3k4me3/Schraham/${samplename}_accepted_hits.bam
    #shortLabel $samplename.bam
    #longLabel $samplename.bam
    #type bam
    #pairEndsByName .
    #pairSearchRange 10000
    #bamColorMode strand
    #maxWindowToDraw 200000
    #db $index
    #visibility pack
    #colorByStrand 200,100,0 0,100,200

    ## ---------- bigwig
    track $samplename.bw
    bigDataUrl http://zlab.umassmed.edu/~dongx/tracks/HD/h3k4me3/Schraham/${samplename}_accepted_hits.bw
    shortLabel $samplename.bw
    longLabel $samplename.bw
    type bigWig
    autoScale on
    maxHeightPixels 100:20:8
    visibility full
    aggregate transparentOverlay
    showSubtrackColorOnUi on
    windowingFunction mean+whiskers
    priority 1.4
    configurable on
    color $color

    " | sed 's/^[ \t]*//g' >> $HOME/projects/bu_neuro/data/trackDb.ChIPseq.txt
done

scp $HOME/projects/bu_neuro/data/trackDb.ChIPseq.txt zlab:~/public_html/tracks/HD




# old version
#controlIP=~/nearline/BU/h3k4me3/2037in.dt1.bed
#macs14 -t Rick_h3k4me3_combined.accepted_hits.bed -c $controlIP -f BED -n Rick_h3k4me3_combined -w -S --call-subpeaks 2> Rick_h3k4me3_combined_macs.log &
#macs14 -t JF_h3k4me3_combined.accepted_hits.bed -c $controlIP -f BED -n JF_h3k4me3_combined -w -S 2> JF_h3k4me3_combined_macs.log &
#macs14 -t Schraham_h3k4me3_control_combined.accepted_hits.bed -c $controlIP -f BED -n Schraham_h3k4me3_control_combined -w -S 2> Schraham_h3k4me3_control_combined_macs.log &
#mkdir vs.2037in; mv Rick_h3k4me3_combined* JF_h3k4me3_combined* vs.2037in/

#####  differential binding between samples
#mkdir RK.vs.Schraham_3706r
#cd !$
#export MAnorm=$HOME/bin/MAnorm_Linux_R_Package/
# $MAnorm/MAnorm2.sh ../vs.2037in/Rick_h3k4me3_combined_peaks.bed ../Schraham_3706r/Schraham_3706r_peaks.bed ../Rick_h3k4me3_combined.accepted_hits.bed ../Schraham_3706r/accepted_hits.bed 200 200 &
#mkdir  ~/scratch/RK.vs.Schrahm_combined;
#cd ~/scratch/RK.vs.Schrahm_combined
#export MAnorm=$HOME/bin/MAnorm_Linux_R_Package/
# $MAnorm/MAnorm2.sh ../vs.2037in/Rick_h3k4me3_combined_peaks.bed ../vs.2037in/Schraham_h3k4me3_control_combined_peaks.bed ../Rick_h3k4me3_combined.accepted_hits.bed ../Schraham_h3k4me3_control_combined.accepted_hits.bed 200 200 &
#
#
#controlIP=~/scratch/Schraham_3706r/accepted_hits.bed
#macs14 -t Rick_h3k4me3_combined.accepted_hits.bed -c $controlIP -f BED -n Rick_h3k4me3_combined -w -S 2> Rick_h3k4me3_combined_macs.log &
#macs14 -t JF_h3k4me3_combined.accepted_hits.bed -c $controlIP -f BED -n JF_h3k4me3_combined -w -S 2> JF_h3k4me3_combined_macs.log &
#mkdir vs.3706r; mv Rick_h3k4me3_combined* JF_h3k4me3_combined* vs.3706r/
#
#controlIP=Schraham_h3k4me3_control_combined.accepted_hits.bed
#macs14 -t Rick_h3k4me3_combined.accepted_hits.bed -c $controlIP -f BED -n Rick_h3k4me3_combined -w -S 2> Rick_h3k4me3_combined_macs.log &
#macs14 -t JF_h3k4me3_combined.accepted_hits.bed -c $controlIP -f BED -n JF_h3k4me3_combined -w -S 2> JF_h3k4me3_combined_macs.log &
#mkdir vs.combinedCt; mv Rick_h3k4me3_combined* JF_h3k4me3_combined* vs.combinedCt/
