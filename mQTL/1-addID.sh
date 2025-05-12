#for chr in $(seq 1 22); do
chr=AFR
cut -f7- -d" " ../Imputed.convert.chr$chr.raw > a.txt

sed -i "s/ /\t/g" a.txt
### need to check the number of colname in colname.txt and in Imputed.convert.chr$chr.raw files
cut -f2 -d" " ../Imputed.convert.chr$chr.raw > a.log
tail -n +2 a.log > b.log
cut -f2 -d"_" b.log > a.log
awk -f ~/KoborLab/kobor_space/zdong/Monkey/Probes3/Old_Genotype-Methylation/Brain_GSE112525/A1-b37.Ilmn_strand_flip/Correlation/transpose1.awk a.log > Genotye_samplenames.txt
rm a.log b.log
## In this dataset, the colname file has been modified to match the nimber

awk -f ~/KoborLab/kobor_space/zdong/Monkey/Probes3/Old_Genotype-Methylation/Brain_GSE112525/A1-b37.Ilmn_strand_flip/Correlation/transpose1.awk a.txt > chr$chr.transpose.txt

#cat colname.txt b.txt > chr$chr.transpose.txt

rm -f a.txt
# done
