i=AFR
plink --noweb --file $i/GSE24274/Imputation/Combined --extract list.snps --make-bed --out data1
plink --noweb --file $i/HapMap3v3/Imputation/Combined --extract list.snps --make-bed --out data2
plink --noweb --file $i/WGBS_impute/Combined --extract list.snps --make-bed --out data3
plink --noweb --merge-list allfiles_$i.txt --make-bed --out merge.$i


## convert ped format into the 012 format (o hom ancesteral), 1 het, 2 dom derived
# for chr in $(seq 1 22); do
plink --bfile merge.$i  --recodeA --out Imputed.convert.chr$i
# done

i=EUR
plink --noweb --file $i/GSE24260/Impute/Combined --extract list.snps --make-bed --out data4
plink --noweb --file $i/HapMap3v3/Imputation/Combined --extract list.snps --make-bed --out data5
plink --noweb --file $i/WGBS_impute/Combined --extract list.snps --make-bed --out data6
plink --noweb --merge-list allfiles_$i.txt --make-bed --out merge.$i
# done

## convert ped format into the 012 format (o hom ancesteral), 1 het, 2 dom derived
# for chr in $(seq 1 22); do
plink --bfile merge.$i  --recodeA --out Imputed.convert.chr$i
# done

## Combine AFR and EUR
plink --noweb --merge-list allfiles_both.txt --make-bed --out merge.both
## convert ped format into the 012 format (o hom ancesteral), 1 het, 2 dom derived
# for chr in $(seq 1 22); do
plink --bfile merge.both  --recodeA --out Imputed.convert.chrboth
plink --bfile merge.both --recode --tab --out merge.both