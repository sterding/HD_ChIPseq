#!/bin/sh
#$ -V
#$ -pe single 8
#$ -o $HOME/sge_jobs_output/sge_job.$JOB_ID.out -j y
#$ -M sterding.hpcc@gmail.com
#$ -S /bin/bash
#$ -m e

ALIGNER=$1

inputdir="/home/dongx/scratch"
BOWTIE_INDEXES="/home/dongx/RDRIVE/hg19"

## --------------------------------
##  Cuffdiff  (only works with tophat output because of the XS:+/- tag)
## --------------------------------
if [[ "$ALIGNER" == "tophat" ]]; then

    [ -d $HOME/scratch/HD_cuffdiff_$ALIGNER ] || mkdir -p $HOME/scratch/HD_cuffdiff_$ALIGNER

    gtflist=""
    samlist=""
    for inputfile1 in `ls -d $inputdir/HD_*{MyersLab,JiangFanChen}*fastq`
    do
        fqfile=$( basename $inputfile1 )
        mv ${inputfile1}/cufflinks_allknown_$ALIGNER/transcripts.gtf $HOME/scratch/HD_cuffdiff_$ALIGNER/$fqfile.cufflinks_allknown_transcripts.gtf
        gtflist="$gtflist $HOME/scratch/HD_cuffdiff_$ALIGNER/$fqfile.cufflinks_allknown_transcripts.gtf"
        fqfile=${fqfile##HD_}  # ${string##substring}: Deletes longest match of $substring from front of $string.
        samlist="$samlist ${inputfile1}/$ALIGNER.$fqfile.sam"
    done

    # to void the path bug in cuffcompare
    cd $HOME/scratch/HD_cuffdiff_$ALIGNER

    cuffcompare -o HD_cuffcompare_$ALIGNER -r $BOWTIE_INDEXES/hg19.protein-coding.gtf $gtflist
    cuffdiff -p 8 -v -L Ct_B6341,HD_B4242,HD_B3584,Ct_B3732,Ct_B6096,HD_B4430,HD_B4189,Ct_B6124,Ct_B6356,HD_B4412 -o ./ HD_cuffcompare_$ALIGNER.combined.gtf $samlist

    # post processing
    mv $HOME/scratch/HD_cuffdiff_$ALIGNER /home/dongx/projects/bu_neuro/results
fi

## --------------------------------
##  htseq
## --------------------------------
# after htseq
[ -d $HOME/scratch/HD_htseq_$ALIGNER ] || mkdir -p $HOME/scratch/HD_htseq_$ALIGNER

echo -e "ID\tCt_B6341\tHD_B4242\tHD_B3584\tCt_B3732\tCt_B6096\tHD_B4430\tHD_B4189\tCt_B6124\tCt_B6356\tHD_B4412"  > $HOME/scratch/HD_htseq_${ALIGNER}/$ALIGNER.hgseqcount.tab
paste ~/scratch/HD_*.fastq/$ALIGNER*hgseqcount*.tab | cut -f1,2,4,6,8,10,12,14,16,18,20 >> $HOME/scratch/HD_htseq_${ALIGNER}/$ALIGNER.hgseqcount.tab

# post processing
cp -R $HOME/scratch/HD_htseq_${ALIGNER} /home/dongx/projects/bu_neuro/results
