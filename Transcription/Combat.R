setwd("~/KoborLab/kobor_space/zdong/Population_dataset/MHB_May3/WGBS_transcript")
library(sva)
test<-read.table("Combined_needadjust_batch.txt",header=T,row.names = 1)
array<-read.table("../Array_transcript/Array_needadjust_batch.txt",header=T,row.names = 1)
a<-read.table("README.txt",header = T)
a$Tans_ID<-paste(a$Tans_ID,".CEL.gz", sep = "")
test1<-test[,c("GSM1202821_b1_lcl_hg133plus2.cel.gz","GSM1202824_b2_lcl_hg133plus2.cel.gz","GSM1202827_b3_lcl_hg133plus2.cel.gz",
               "GSM291698.CEL.gz","GSM291618.CEL.gz",
               as.character(a$Tans_ID))]
#colnames(test1)<-c("SRR948847","SRR948848","SRR948849","SRR568016","SRR4453293",
#                   "SRR1282079","SRR1282090","SRR1282094","SRR1282097","SRR1282101","SRR1282105","SRR1282109","SRR1282111","SRR1867805")
test1<-cbind(test1,array)
subgroups<-c(replicate(60,1),replicate(54,0),replicate(96,1),replicate(95,0))
batch<-c(1,1,1,replicate(111,2),replicate(191,0))
pheno<-data.frame(subgroups,batch)
rownames(pheno)<-colnames(test1)
write.table(pheno,"pheno.txt",sep="\t",quote = F)
modcombat<-model.matrix(~1, data=pheno)
combat_mydata= ComBat(dat=as.matrix(test1), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
total<-data.frame(test[,1:2],combat_mydata)
write.table(total,file="combat_total_WGBSandArray.txt",quote=F,sep="\t",col.names = T,row.names = T)
# total1<-total[,1:16]
# write.table(total1,file="combat_total_witharray.txt",quote=F,sep="\t",col.names = T,row.names = T)
# total2<-total[,c(1:2,17:ncol(total))]
# write.table(total2,file="../Array_transcript/combat_array_withWGBS.txt",quote=F,sep="\t",col.names = T,row.names = T)

