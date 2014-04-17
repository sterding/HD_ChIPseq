Analysis pipeline for the HD H3K4me3 project
============================================

# 1. reads mapping and peak calling for all H3K4me3 ChIP-seq samples, incl. healthy control and HD
make_sge_ChIPseq.sh

for i in ~/nearline/BU/h3k4me3/h3k4me3_HD_batch1/*.fq; do j=`basename $i`; qsub ~/projects/bu_neuro/src/make_sge_ChIPseq.sh $i h3k4me3_HD_batch1_${j/.*/}; done
for i in ~/nearline/BU/h3k4me3/control_from_Schraham/*.fq; do j=`basename $i`; qsub ~/projects/bu_neuro/src/make_sge_ChIPseq.sh $i Schraham_${j/.*/}; done

# 2. merge peaks called from individual alignment
peakCalling.sh

Two stratedgies:
1) merge alignment and then call peaks from the merged alignment
2) merge peaks called from individual alignment. Only a peak occurred in >=3 samples will be defined as a peak.

It shows that the 2nd stratedgy is more robust. (see figure: https://docs.google.com/presentation/d/1kesRVX4McijTVSfk2ul6c_HrSe88pNuEaVzeRq8B8Dc/edit#slide=id.g12433e2fc_05)

Output of 2) method:
Table 2: https://docs.google.com/presentation/d/1kesRVX4McijTVSfk2ul6c_HrSe88pNuEaVzeRq8B8Dc/edit#slide=id.g12433e2fc_049

# 3. top 1000 H3k4me3 peaks in HD and healthy controls, dividing into promoter-proximal and promoter-distal peaks for each group

# 4. genes with differential H3k4me3 in proximal promoter regions
DP_ChIPseq.sh

# 4.1 pathway analysis
# 4.2 line up with gene with differential expression from RNAseq

# 5. enhancers with differential H3k4me3 signal
DE_ChIPseq.sh

# 5.1 the nearest genes and
# 5.2 its pathway analysis