cd /Users/xdong/projects/HD/data/scratch/

enhancers_CAGE=~/projects/HD/data/hg19_permissive_enhancers_expression_rle_tpm.brain.tab  # in total, 43011 enhancers defined by CAGE
# wget http://enhancer.binf.ku.dk/presets/hg19_permissive_enhancers_expression_rle_tpm.csv.gz
# rowsToCols hg19_permissive_enhancers_expression_rle_tpm.csv -fs=, stdout | grep -E "chr1:858256-858648|CNhs10617|CNhs10647" | rowsToCols -tab stdin stdout | sed 's/"//g' | awk '{OFS="\t"; print $3, $1,$2;}' | sed 's/[:-]/ \t/g' | grep chr > hg19_permissive_enhancers_expression_rle_tpm.brain.tab
# add header "#chr	start	end	TPM_frontallobeadult_CNhs10647	TPM_brainadult_CNhs10617"

# get max normalized signal of H3k4me3 per brain enhancers (defined by CAGE)
awk '{OFS="\t"; if($4>0) print $1,$2,$3,$1"_"$2"_"$3"_"$4;}' $enhancers_CAGE | bigWigAverageOverBed  trimmedmean.normalized.Ct.bw stdin stdout | cut -f1,4 | sed 's/_/\t/g' > enhancers_CAGE.h3k4me3inHC.tab

awk '$4>0 || $5>0' $enhancers_CAGE > test.bed  # 12097 CAGE enhancers expressed in brain/frontal lobe

~/projects/HD/src/toBinRegionsOnBigwig.sh trimmedmean.normalized.Ct.bw test.bed 1 | paste - test.bed > enhancers_CAGE.expressedinBrain.h3k4me3inHC.tab 

grep distal h3k4me3_peaks_intersected.Ct.bed.commonuniq.inpromoter | grep common | intersectBed -b - -a test.bed -u > test.bed2
~/projects/HD/src/toBinRegionsOnBigwig.sh trimmedmean.normalized.Ct.bw test.bed2 1 | paste - test.bed2 | sort -k2,2nr > enhancers_CAGE.HCdistalcommon.tab

grep distal h3k4me3_peaks_intersected.HD.bed.commonuniq.inpromoter | grep uniq | awk '($3-$2)>100' | intersectBed -b - -a test.bed -u > test.bed2
~/projects/HD/src/toBinRegionsOnBigwig.sh trimmedmean.normalized.HD.bw test.bed2 1 | paste - test.bed2 | sort -k2,2nr > enhancers_CAGE.HDdistaluniq.tab 

grep distal h3k4me3_peaks_intersected.Ct.bed.commonuniq.inpromoter | grep uniq | awk '($3-$2)>100' | intersectBed -b - -a test.bed -u > test.bed2
~/projects/HD/src/toBinRegionsOnBigwig.sh trimmedmean.normalized.Ct.bw test.bed2 1 | paste - test.bed2 | sort -k2,2nr > enhancers_CAGE.HCdistaluniq.tab 

R
df=read.table("enhancers_CAGE.HDdistaluniq.tab")
cor(df$V2,df$V6, method="spearman")
[1] 0.3404799

df[df$V2>1 & df$V6>1,]

df=read.table("data/hg19_permissive_enhancers_expression_rle_tpm.brain.tab.overlapwithh3k4me3distalpeakinHC.tab")
cor(df$V2,df$V6, method="spearman")
[1] 0.3404799

## union of CAGE enhancer and H3k4me3 distal peaks
##------------------------------------------------
cat *commonuniq.inpromoter | grep distal | bedSort stdin stdout | mergeBed | intersectBed -a - -b $enhancers_CAGE -u | wc -l
# 1208
cat *commonuniq.inpromoter | grep "common:distal" | bedSort stdin stdout | mergeBed | intersectBed -a - -b $enhancers_CAGE -u | wc -l
# 957

cat *commonuniq.inpromoter | grep distal | intersectBed -a $enhancers_CAGE -b - -wo | awk '{OFS="\t"; print $1, $2,$3; print $6,$7,$8;}' | bedSort stdin stdout | mergeBed > enhancers_CAGE.plus.distalpeaks.bed
# 1095

paste <(~/projects/HD/src/toBinRegionsOnBigwig.sh trimmedmean.normalized.Ct.bw enhancers_CAGE.plus.distalpeaks.bed 1) <(~/projects/HD/src/toBinRegionsOnBigwig.sh trimmedmean.normalized.HD.bw enhancers_CAGE.plus.distalpeaks.bed 1 | cut -f2 -d' ')  > enhancers_CAGE.plus.distalpeaks.bed.h3k4me3summit.in.HCHD

#R
df = read.table("enhancers_CAGE.plus.distalpeaks.bed.h3k4me3summit.in.HCHD")
colnames(df)=c("coordinate","HC","HD")
df=cbind(df, M=log2(df$HD/df$HC), A=(df$HD+df$HC)/2)
pdf("../../result/2014Jun12/enhancers_CAGE.plus.distalpeaks.bed.h3k4me3summit.MAplot.pdf", width=6, height=6)
plot(M~A, df, pch=21, bg=ifelse(abs(df$M)>=1,'red','#222222'), col='white', cex=.7, xlab="Mean of submit of H3K4me3 normalized signal", ylab="foler change of HD vc. HC (log2)")
abline(h=0, lty=2, lwd=1)
dev.off()


