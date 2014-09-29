Analysis pipeline for the HD H3K4me3 project
============================================

# 1. reads mapping and peak calling for all H3K4me3 ChIP-seq samples, incl. healthy control and HD
make_sge_ChIPseq.sh

for i in ~/nearline/BU/h3k4me3/h3k4me3_HD_batch1/*.fq; do j=`basename $i`; qsub ~/projects/bu_neuro/src/make_sge_ChIPseq.sh $i h3k4me3_HD_batch1_${j/.*/}; done
for i in ~/nearline/BU/h3k4me3/control_from_Schraham/*.fq; do j=`basename $i`; qsub ~/projects/bu_neuro/src/make_sge_ChIPseq.sh $i Schraham_${j/.*/}; done

* Note: Since we are not going to use batch2 and JF's data, we won't talk about it here.

# 2. merge peaks called from individual alignment
peakCalling.sh

Two stratedgies including:
1) merge alignment and then call peaks from the merged alignment
2) merge peaks called from individual alignment. Only a peak occurred in >=3 samples in each category will be defined as a peak. 

and it shows that the 2nd stratedgy is more robust. (see figure: https://docs.google.com/presentation/d/1kesRVX4McijTVSfk2ul6c_HrSe88pNuEaVzeRq8B8Dc/edit#slide=id.g12433e2fc_05)

Output of 2) method:
Table 2: https://docs.google.com/presentation/d/1kesRVX4McijTVSfk2ul6c_HrSe88pNuEaVzeRq8B8Dc/edit#slide=id.g12433e2fc_049

$ wc -l *bed
28579 h3k4me3_peaks_intersected.HD.bed
33144 h3k4me3_peaks_intersected_Ct.bed

# 3. revised method:
peakAnalysis.sh

1) take union of peaks from HD and control and merge them into one big set of peaks;
2) for each peak in the set, measure its distance from middle point to the nearest TSS in Gencode annotation (v17), normalized RPM signal in HD and control, # of samples supported that peak. For peak with distance less or equal to 500bp from TSS, it's annotated as a proximal peak; otherwise, distal one. 

output file: h3k4me3_peaks_union.signal.neighborhood.annotated.xls (N=28608)
  
# 3. significant H3k4me3 peaks in HD and healthy controls, dividing into 2x2 categroies (uniq/common and proximal/distal)
peakAnalysis.sh

Output: a list of differental/uniq peaks and their nearest genes, 2-3 cases study

# 4. genes with differential H3k4me3 in proximal promoter regions
# use reads count in [-500,+500] of TSS to call DESeq2
DP_ChIPseq.sh

## 4.1 pathway analysis
## 4.2 line up with gene with differential expression from RNAseq

# 5. brain-specific enhancers with differential H3k4me3 signal
DE_ChIPseq.sh

## 5.1 the nearest genes and possible external support of their association (e.g. ChIA-PET, DNaseI, latest FANTOM5 enhancer-TSS association data etc.)
## 5.2 GO / pathway analysis for the nearest genes