#setwd("~/KoborLab/kobor_space/kandy/home/zdong/Population_dataset/MHBmap/SNP_new_forimputation/Correation/AFR")
##================================Generate methylation input file======================================
meth<-read.table("../../Array_methylation/MHB_recall/AML/Call_Pop-MHBs_in_WGBS/AML_arrayinWGBS_wilcox.txt",header = T,row.names = 1)
meth_names<-read.table("../../WGBS_methylation/sampleInfo.txt")
#sort -k2,2 -k3,3n chrAFR.transpose.txt > all.transpose_sorted.txt
#remove all colnames 
# cat Genotye_samplenames.txt all.transpose_sorted.txt > all.transpose.txt
genotype<-read.table("all.transpose.txt",header=T,row.names = 1)
### get the list from SNP header
genotype_names<-meth_names
rownames(genotype_names)<-sub("[.].*", "", rownames(genotype_names))
rownames(genotype_names)<-sub("fileID_", "", rownames(genotype_names))
genotype_names<-genotype_names[1:10,]

share<-intersect(meth_names$Sample,genotype_names$Sample)
meth_names<-meth_names[meth_names$Sample %in% share,]
meth_names<-meth_names[order(meth_names$Sample),]
genotype_names<-genotype_names[genotype_names$Sample %in% share,]
genotype_names<-genotype_names[order(genotype_names$Sample),]
meth1<-meth[,rownames(meth_names)]
colnames(meth1)<-meth_names$Sample
write.table(meth1,file = "Input_methylation.txt",quote=F,sep="\t")
# add "TargetID"

genotype1<-genotype[,rownames(genotype_names)]
colnames(genotype1)<-genotype_names$Sample
write.table(genotype1,file = "Input_genotype_allchr.txt",quote=F,sep="\t")
# add "TargetID"

##================================Generate covariates input file======================================
confounders<-meth_names
colnames(confounders)[1]<-"TargetID"
confounders1<-t(confounders)
confounder2<-confounders1[1:3,] ##!!!!! only keep ethnicity into our analysis
write.table(confounder2,file="Input_covariate.txt",quote = F,sep="\t",col.names = F)




