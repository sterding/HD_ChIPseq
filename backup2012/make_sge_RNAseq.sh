#!/bin/sh

# Usage: ./make_sge_RNAseq.sh $HOME/nearline/BU/Jul2012/rnaseq/filtered

cpu=24
index=hg19
inputdir=$1  #$HOME/nearline/BU/Jul2012/rnaseq/filtered
outputdir=$HOME/scratch
adaptorfile=adaptor.fa
ANNOTATION=$GENOME/hg19/Annotation/Genes
BOWTIE_INDEXES=$GENOME/hg19/Sequence/BowtieIndex

###########################################
#============= mapping options
###########################################

#phred
bowtie="--phred33-quals"; bowtie2="--phred33"; tophat=""; far="fastq-sanger"; fastqmcf="33"; trimmomatic="-phred33"
#mismatch
mm=2
#PE option
PE="--mate-inner-dist 50"
#strand
strandoption="--library-type fr-unstranded"


cd $inputdir
[ -e $HOME/projects/bu_neuro/data/trackDb.txt ] && rm $HOME/projects/bu_neuro/data/trackDb.txt;

c=0;h=0;gtflist="";samlist=""; labels=""

for i in Sample_[!u]*.R1.fastq.gz;
do
    R1=$i
    R2=${i/R1/R2};
    samplename=${R1/.R*/}

    #echo $R1, $R2

    if [[ "$samplename" =~ "_H_" ]]; then
        let h=$h+20
        color="255,"$h",51"
    else
        let c=$c+20
        color=$c",255,0"
    fi

    ### generate trackDb.txt
    #echo "
    #######################################################
    ###  $samplename
    #######################################################
    ### ---------- BAM
    ##track $samplename.bam
    ##bigDataUrl http://zlab.umassmed.edu/~dongx/tracks/HD/$samplename.accepted_hits.bam
    ##shortLabel $samplename.bam
    ##longLabel $samplename.bam
    ##type bam
    ##pairEndsByName .
    ##pairSearchRange 10000
    ##bamColorMode strand
    ##maxWindowToDraw 200000
    ##db $index
    ##visibility pack
    ##colorByStrand 200,100,0 0,100,200
    #
    ### ---------- bigwig
    #track $R1.bigwig
    #bigDataUrl http://zlab.umassmed.edu/~dongx/tracks/HD/$samplename.accepted_hits.bw
    #shortLabel $samplename.bigwig
    #longLabel $samplename.bigwig
    #type bigWig
    #autoScale on
    #visibility full
    #aggregate transparentOverlay
    #showSubtrackColorOnUi on
    #windowingFunction mean+whiskers
    #priority 1.4
    #configurable on
    #color $color
    #
    #" | sed 's/^[ \t]*//g' >> $HOME/projects/bu_neuro/data/trackDb.txt

    [ -d $outputdir/$samplename ] || mkdir -p $outputdir/$samplename

    echo "
    ###########################################
    ## bash script for running RNAseq QC on cluster
    ###########################################
    #!/bin/sh
    #$ -V
    #$ -pe openmpi $cpu
    #$ -cwd
    #$ -o \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y
    #$ -S /bin/bash
    #$ -l mem_free=6G

    ###########################################
    ############## 1. Configuring
    ###########################################

    export BOWTIE_INDEXES=\$GENOME/$index/Sequence/Bowtie2Index/
    export BOWTIE2_INDEXES=\$GENOME/$index/Sequence/Bowtie2Index/
    export ANNOTATION=\$GENOME/$index/Annotation/Genes

    ln -sf \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out $outputdir/$samplename/sge.log

    cd $inputdir

    ###########################################
    ###############  2. quality filter: adaptor removal/clip
    ###########################################

    ##### adaptor removal
    [ -d filtered ] || mkdir filtered
    ##java -classpath \$CLASSPATH/trimmomatic-0.22.jar org.usadellab.trimmomatic.TrimmomaticPE -threads $cpu $trimmomatic $R1 $R2 filtered/$R1 filtered/${R1/R1/R1.unpaired} filtered/$R2 filtered/${R2/R2/R2.unpaired} LEADING:3 TRAILING:3 ILLUMINACLIP:$adaptorfile:2:40:15 SLIDINGWINDOW:4:15 MINLEN:16;
    fastq-mcf -o filtered/$R1 -o filtered/$R2 -l 16 -q 15 -w 4 -x 10 -u -P $fastqmcf $adaptorfile $R1 $R2

    cd filtered

    #############################################
    ################ 3. QC
    ############################################

    fastqc --outdir $outputdir/$samplename --extract -t 2 $R1 $R2
    rm $outputdir/$samplename/*fastqc.zip

    ############################################
    ################ 4. mapping to the genome
    ############################################
    #
    ## tophat (output accepted_hits.sam, allow multiple hits)
    tophat -o $outputdir/$samplename --no-convert-bam -p $cpu --read-mismatches $mm $tophat $PE $strandoption --min-anchor-length 8 --min-intron-length 30 --max-intron-length 50000 --splice-mismatches 1 --max-multihits 100 --no-coverage-search genome_offrate3 $R1 $R2

    ###########################################
    ############### 5. post-processing, format converting
    ###########################################


    cd $outputdir/$samplename

    ## sam -> bam -> sorted -> index
    samtools view -Sbut \$BOWTIE_INDEXES/genome.fai accepted_hits.sam | samtools sort - accepted_hits.sorted
    mv accepted_hits.sorted.bam accepted_hits.bam
    samtools index accepted_hits.bam


    ###########################################
    ################# 7. assembly
    ###########################################

    cd $outputdir/$samplename

    echo \"## run cufflinks to get FPKM\"
    cufflinks -q --no-update-check $strandoption -o ./ -p $cpu -G \$ANNOTATION/gencode.v14.annotation.gtf -M \$ANNOTATION/chrM.rRNA.tRNA.gtf -b \$BOWTIE_INDEXES/genome.fa --multi-read-correct accepted_hits.bam
    #cufflinks -q --no-update-check $strandoption -o ./ -p $cpu -G \$ANNOTATION/genes.gtf -M \$ANNOTATION/chrM.rRNA.tRNA.gtf -b \$BOWTIE_INDEXES/genome.fa --multi-read-correct accepted_hits.bam

    #echo \"## run htseq for reads count\"
    #htseq-count -m intersection-strict -t exon -i gene_id -s no -q accepted_hits.sam \$ANNOTATION/gencode.v14.annotation.gtf > $samplename.hgseqcount.by.gene.tab 2> $samplename.hgseqcount.by.gene.tab.stderr
    #echo \"## run bedtools for reads count\"
    #bedtools multicov -D -split -bams accepted_hits.bam -bed \$ANNOTATION/gencode.v14.annotation.bed15 > $samplename.bedtools.by.trans.tab

    ############################################
    ############### 6. prepare for tracks files to display on UCSC
    ############################################
    #
    #[ -d $outputdir/ucsc ] || mkdir $outputdir/ucsc
    cd $outputdir/ucsc

    ## make index for the (sorted) BAM
    #mv $outputdir/$samplename/accepted_hits.bam $samplename.accepted_hits.bam
    #mv $outputdir/$samplename/accepted_hits.bam.bai $samplename.accepted_hits.bam.bai
    #
    ## QC
    #mv $outputdir/$samplename/*_fastqc ./

    ## bigwig
    #bamToBed -i $samplename.accepted_hits.bam -split > $samplename.accepted_hits.bed
    #sort -k1,1 $samplename.accepted_hits.bed | bedItemOverlapCount mm9 -chromSize=\$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n > $samplename.accepted_hits.bedGraph
    #bedGraphToBigWig $samplename.accepted_hits.bedGraph \$ANNOTATION/ChromInfo.txt $samplename.accepted_hits.bw
    #rm $samplename.accepted_hits.bed $samplename.accepted_hits.bedGraph

    cp $outputdir/$samplename/isoforms.fpkm_tracking $samplename.isoforms.fpkm_tracking
    cp $outputdir/$samplename/genes.fpkm_tracking $samplename.genes.fpkm_tracking

    # gtf of assembly
    echo \"track name=$samplename description=$samplename visibility=pack colorByStrand='200,100,0 0,100,200'\" > $samplename.transcripts.gtf
    cat $outputdir/$samplename/transcripts.gtf >> $samplename.transcripts.gtf
    gzip -f $samplename.transcripts.gtf

    echo \"JOBDONE!\"

    " | sed 's/^[ \t]*//g'  > $outputdir/$samplename/$samplename.sge
    #cat $outputdir/$samplename/$samplename.sge
    jobid=`qsub $outputdir/$samplename/$samplename.sge | cut -f3 -d' '`

    echo "Your job is submitted (jobID: $jobid) with SGE script at $outputdir/$samplename/$samplename.sge"

    gtflist="$gtflist $outputdir/$samplename/transcripts.gtf"
    samlist="$samlist $outputdir/$samplename/accepted_hits.sam"
    if [ "$labels" == "" ];
    then
        labels="$samplename";
    else
        labels="$labels,$samplename"
    fi

done

exit;


# run cuffdiff

[ -d $outputdir/HD_cuffdiff ] || mkdir -p $outputdir/HD_cuffdiff

echo "
###########################################
## bash script for running cuffdiff on cluster
###########################################
#!/bin/sh
#$ -V
#$ -pe openmpi $cpu
#$ -cwd
#$ -o \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y
#$ -S /bin/bash
#$ -l mem_free=3G

export BOWTIE2_INDEXES=\$GENOME/$index/Sequence/Bowtie2Index/
export ANNOTATION=\$GENOME/$index/Annotation/Genes

[ -d $outputdir/HD_cuffdiff ] || mkdir -p $outputdir/HD_cuffdiff
cd $HOME/scratch/HD_cuffdiff

#cuffcompare
cuffcompare -s \$GENOME/hg19/Sequence/Chromosomes -r $ANNOTATION/genes.gtf $gtflist

# HD vs. Control
HD=`echo $samlist | sed 's/ /\n/g' | grep \"_H_\" | tr '\n' ',' | sed 's/,$//g'`
Ct=`echo $samlist | sed 's/ /\n/g' | grep \"_C_\" | tr '\n' ',' | sed 's/,$//g'`
cuffdiff -p 8 -v -L HD,Ct $strandoption -M \$ANNOTATION/chrM.rRNA.tRNA.gtf -b \$BOWTIE2_INDEXES/genome.fa -u -o ./cuffdiff2 cuffcmp.combined.gtf \$HD \$Ct &
cuffdiff -p 8 -v -L HD,Ct $strandoption -M \$ANNOTATION/chrM.rRNA.tRNA.gtf -b \$BOWTIE2_INDEXES/genome.fa -u -o ./cuffdiff3 \$ANNOTATION/genes.gtf \$HD \$Ct &

# all to all
#cuffdiff -p 8 -v -L $labels $strandoption -M \$ANNOTATION/chrM.rRNA.tRNA.gtf -b \$BOWTIE2_INDEXES/genome.fa -u cuffcmp.combined.gtf $samlist
cuffdiff -p 8 -v -L $labels $strandoption -M \$ANNOTATION/chrM.rRNA.tRNA.gtf -b \$BOWTIE2_INDEXES/genome.fa -u \$ANNOTATION/genes.gtf $samlist


" | sed 's/^[ \t]*//g'  > $outputdir/HD_cuffdiff/submit.sge

cat $outputdir/HD_cuffdiff/submit.sge

jobid=`qsub $outputdir/HD_cuffdiff/submit.sge | cut -f3 -d' '`
ln -fs $HOME/sge_jobs_output/sge_job.$jobid.out $outputdir/HD_cuffdiff/sge.log

echo "
Your job is submitted (jobID: $jobid) !
The SGE script is saved as $outputdir/HD_cuffdiff/submit.sge
The running log is saved as $outputdir/HD_cuffdiff/sge.log
"
