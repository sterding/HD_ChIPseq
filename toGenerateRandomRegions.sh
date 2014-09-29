# ===============================================================
# Script to generate a random background regions with the same length distrubiton as input bed
# Author: Xianjun Dong
# Date: Jun 3, 2014
# Usage:
# ../../src/toGenerateRandomRegions.sh input.bed > output.random.bed
# or
# cat input.bed | ../../src/toGenerateRandomRegions.sh -
# Requirement:
# 1. input bed file at least have chr, start, and end fields
# ===============================================================

bedfile=$1

# download hg.genome
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo"  > hg19.genome

cat $bedfile | while read chr start end rest
do
    let l=$end-$start;
    bedtools random -g hg19.genome -n 1 -l $l
done