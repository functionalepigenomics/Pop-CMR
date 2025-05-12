chr=AFR
cut -f1,4  ../merge.$chr.map > a.log
#cut -f1 chr$chr.transpose.txt | sed '1d' > b.log
cut -f1 chr$chr.transpose.txt > b.log
paste a.log b.log > SNP_coordinate_chr$chr.txt
rm -f a.log b.log
#done
cat SNP_coordinate_chr*.txt > SNP_coordinate.txt
rm SNP_coordinate_chr*
sed -i -e 's/^/chr/g' SNP_coordinate.txt
#cut -f1 ../AFR_datanew.txt > b.log
#awk 'FNR==NR{a[$1];next}($1 in a){print}' b.log ~/KoborLab/kobor_space/kandy/home/zdong/Monkey/Probes2/Mdata/IlluminaHumanMethylation450kanno.ilmn12.hg19.txt | cut -f1-3 > CpG_coordinate.txt
#rm b.log

Rscript SNPcoordination.R

#sed -i -e 's/^/chr/g' AFR_CpG_coordinate1.bed
sort -k2,2 -k3,3n MHB_coordinate1.bed > CpG_coordinate.sorted1.bed
sort -k2,2 -k3,3n SNP_coordinate1.bed > SNP_coordinate.sorted1.bed

#bedtools intersect -a CpG_coordinate.sorted.bed -b SNP_coordinate.sorted.bed -wa -wb > CpG-SNP.txt
#bedmap --echo --echo-map --delim '\t' --skip-unmapped CpG_coordinate.sorted.bed SNP_coordinate.sorted.bed > CpG-SNP.txt

