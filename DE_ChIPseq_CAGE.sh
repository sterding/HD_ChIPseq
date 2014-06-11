enhancers_CAGE=~/projects/HD/data/hg19_permissive_enhancers_expression_rle_tpm.brain.tab  # in total, 43011 enhancers defined by CAGE
# wget http://enhancer.binf.ku.dk/presets/hg19_permissive_enhancers_expression_rle_tpm.csv.gz | 
# rowsToCols hg19_permissive_enhancers_expression_rle_tpm.csv -fs=, stdout | grep -E "chr1:858256-858648|CNhs10617|CNhs10647" | rowsToCols -tab stdin stdout | sed 's/"//g' | awk '{OFS="\t"; print $3, $1,$2;}' | sed 's/[:-]/ \t/g' | grep chr > hg19_permissive_enhancers_expression_rle_tpm.brain.tab
# add header "#chr	start	end	TPM_frontallobeadult_CNhs10647	TPM_brainadult_CNhs10617"

cd /Users/xdong/projects/HD/


# get mean normalized signal per brain enhancers (defined by CAGE)
awk '{OFS="\t"; if($4>0) print $1,$2,$3,$1"_"$2"_"$3"_"$4;}' $enhancers_CAGE | bigWigAverageOverBed  data/scratch/trimmedmean.normalized.Ct.bw stdin stdout | cut -f1,4 | sed 's/_/\t/g' > $enhancers_CAGE.h3k4me3inHC.tab

awk '$4>0' $enhancers_CAGE > test.bed

src/toBinRegionsOnBigwig.sh data/scratch/trimmedmean.normalized.Ct.bw test.bed 1 | paste - test.bed > $enhancers_CAGE.h3k4me3inHC.tab 

R
df=read.table("data/hg19_permissive_enhancers_expression_rle_tpm.brain.tab.h3k4me3inHC.tab")
cor(df$V2,df$V6, method="spearman")
[1] 0.3356477

