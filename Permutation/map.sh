#sort-bed /users/zhengd/scratch/Work2/Population/MHBmap/MHB.ethnicity/MHB.GWBSandRRBS/unique.sorted.nonoverlap.YRI.MHB.bed
sort-bed ../Array-results/dmcmr_pos.bed > a.txt
sort-bed ~/KoborLab/kobor_space/zdong/Population_dataset/MHB_New/Analysis/twoloci_cor/Permutation/Permutation/cmr.bed > b.txt
files=(*.bed)
for file in ${files[@]};
do
echo $file
sort-bed $file > sorted.$file.txt
bedmap --echo --echo-map --delim '\t' --skip-unmapped a.txt sorted.$file.txt > $file.cmr_pair.txt

bedmap --echo --echo-map --delim '\t' --skip-unmapped a.txt sorted.$file.txt | wc -l >> cmr_pair.count.txt
for ((i=1;i<=1000000;i++));
do
shuf -n 101 b.txt > n.log
sort-bed n.log > $i.log
bedmap --echo --echo-map --delim '\t' --skip-unmapped $i.log sorted.$file.txt | wc -l >> $file.cmr_background.count.txt
done
rm *.log
done

