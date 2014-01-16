#!/bin/sh
#$ -V
#$ -pe single 8
#$ -o $HOME/sge_jobs_output/sge_job.$JOB_ID.out -j y
#$ -M sterding.hpcc@gmail.com
#$ -S /bin/bash
#$ -m e

inputfile=$1
ALIGNER=$2

# for debug purpose
#inputfile="18_1_62LPFAAXX_092810_MyersLab_B6341.fastq"
#ALIGNER='bwa'

inputdir="/home/dongx/nearline/rnaseq/BU/fastq"
BOWTIE_INDEXES="/home/dongx/RDRIVE/hg19"  # index for both BOWTIE and BWA
OUTPUT_DIR="$HOME/scratch/"

# Make the directory for the job ID you are running if it does not exist
[ -d $HOME/scratch/HD_$inputfile ] || mkdir -p $HOME/scratch/HD_$inputfile

cd $HOME/scratch/HD_$inputfile

#######################################
########### ALIGNMENT #################
#######################################

# 1. Run tophat with strand info  # allow 2 mismatch, unique
if [[ "$ALIGNER" == "tophat" ]]; then
    tophat -p 8 --max-multihits 2 --segment-mismatches 2 --library-type fr-unstranded -o ./ $BOWTIE_INDEXES/hg19 $inputdir/$inputfile
    mv accepted_hits.bam tophat.$inputfile.bam
    mv junctions.bed tophat.$inputfile.junctions.bed

    ## BAM -> SAM:: use -q 20 to filter reads with 0.01 probability of being wrongly aligned.
    samtools view -q 20 tophat.$inputfile.bam -o tophat.$inputfile.sam

    ## only for temporary use
    #cp ~/scratch/oldHD/HD_$inputfile/tophat.mm2.unique.hg19.bam tophat.$inputfile.bam
    #cp ~/scratch/oldHD/HD_$inputfile/tophat.mm2.unique.hg19.sam tophat.$inputfile.sam
    #cp ~/scratch/oldHD/HD_$inputfile/tophat.mm2.unique.hg19.junctions.bed tophat.$inputfile.junctions.bed

    # ----------- generate BAM track in zhome  -----------
    #samtools sort tophat.$inputfile.bam tophat.$inputfile.sorted
    #samtools index tophat.$inputfile.sorted.bam

    # ----------- run cufflinks to get FPKM -----------  CUFFLINKS only works with sam from TOPHAT output!!!!!!!
    # without annotation, used for noval transcription discovery
    #cufflinks -o ./cufflinks -p 8 -M $BOWTIE_INDEXES/hg19.rRNA_MtRNA.gtf -r $BOWTIE_INDEXES/hg19.fa tophat.$inputfile.sam
    # run cuffcompare to compare with reference.gtf
    #cuffcompare -o cuffcompare -r $BOWTIE_INDEXES/hg19.gtf -R -s $BOWTIE_INDEXES/ucsc.chr ./cufflinks/transcripts.gtf

    #  only on hg19.gtf protein-coding genes
    cufflinks -o cufflinks_allknown_$ALIGNER -p 8 -G $BOWTIE_INDEXES/hg19.protein-coding.gtf -r $BOWTIE_INDEXES/hg19.fa -q tophat.$inputfile.sam
    cuffcompare -o cuffcompare_allknown_$ALIGNER -r $BOWTIE_INDEXES/hg19.gtf -R -s $BOWTIE_INDEXES/ucsc_chr ./cufflinks_allknown_$ALIGNER/transcripts.gtf

    # only on lincRNA
    #cufflinks -o ./cufflinks_lincRNA -j 0.05 -p 8 -Q 0 -G $BOWTIE_INDEXES/hg19.lincRNA.gtf -M $BOWTIE_INDEXES/hg19.rRNA_MtRNA.gtf -r $BOWTIE_INDEXES/hg19.fa tophat.$inputfile.sam

fi

# 2. Run BWA with MAPQ<20 filter
if [[ "$ALIGNER" == "bwa" ]]; then
    bwa aln -q 20 $BOWTIE_INDEXES/hg19.fa $inputdir/$inputfile > bwa.$inputfile.sai 2> bwa.$inputfile.sai.stderr
    bwa samse $BOWTIE_INDEXES/hg19.fa bwa.$inputfile.sai $inputdir/$inputfile > bwa.$inputfile.sam 2> bwa.$inputfile.samse.stderr

    # ----------- generate BAM track in zhome  -----------
    #samtools view -bS -o bwa.$inputfile.bam bwa.$inputfile.sam
    #samtools sort bwa.$inputfile.bam bwa.$inputfile.sorted
    #samtools index bwa.$inputfile.sorted.bam
fi

# ----------- move tracks to zhome  -----------
#[ -d $HOME/scratch/HD_tracks ] || mkdir -p $HOME/scratch/HD_tracks
#mv *.$inputfile.sorted.bam $HOME/scratch/HD_tracks/
#mv *.$inputfile.sorted.bam.bai $HOME/scratch/HD_tracks/
# TODO: scp *.bam, *.bam.bai to zhome tracks folder
# rsync -azv $HOME/scratch/HD_tracks dongx@zhome.umassmed.edu:~/public_html


#######################################
########### Assembly  #################
#######################################
# to get reads count on each feature (e.g.gene, transcript). Below use gene
htseq-count -m intersection-strict -t exon -i gene_id -s no -q $ALIGNER.$inputfile.sam $BOWTIE_INDEXES/hg19.protein-coding.gff > $ALIGNER.$inputfile.hgseqcount.by.gene.tab 2> $ALIGNER.$inputfile.hgseqcount.by.gene.tab.stderr
