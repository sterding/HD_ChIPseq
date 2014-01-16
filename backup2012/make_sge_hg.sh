#!/bin/sh

inputdir="/home/dongx/nearline/rnaseq/BU/fastq"

cd $inputdir
for inputfile in `dir -d *{MyersLab,JiangFanChen}*.fastq`
do
    qsub align.sh $inputfile tophat
    qsub align.sh $inputfile bwa
done
