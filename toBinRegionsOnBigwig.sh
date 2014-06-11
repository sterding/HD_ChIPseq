# ===============================================================
# Script for binning a list of regions into N bins and get max bigwig signal on the bins
# Author: Xianjun Dong
# Date: Jun 3, 2014
# Usage:
# toBinRegionsOnBigwig.sh signal.bw input.bed 100
# or
# toBinRegionsOnBigwig.sh http://path.com/signal.bw input.bed 100
# Requirement:
# 1. input bed file at least have chr, start, and end fields
# 2. bigWigSummary (from Jim Kent) should be in the $PATH
# ===============================================================

bigwig=$1
bedfile=$2
N=$3

while read chr start end rest
do
    s=`bigWigSummary $bigwig -udcDir=/tmp -type=max $chr $start $end $N 2>/dev/null | sed 's/n\/a/0/g'`;
    [[ $s == "" ]] && s=`yes 0 | head -n$N | tr '\n' '\t'`
    echo "$chr:$start-$end" $s;
done < $bedfile
