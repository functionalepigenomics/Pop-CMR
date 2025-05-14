data<-read.table(gzfile("beta.txt"),header=T,row.names = 1,nrows = 1000000)
names(data)<-sub("[.]", "", names(data))
names(data)<-sub("[.]", "", names(data))
names(data)<-sub(".CpG_report.merged_CpG_evidence.cov", "", names(data))
data<-data[,!(names(data) %in% c("SRR1282094","SRR1282101","SRR1282105","SRR17509884","SRR948855"))]
ethnicity<-read.table("ethnicity")
afr<-data[,ethnicity$V1=="AFR"] # 54
eur<-data[,ethnicity$V1=="EUR"] # 8
p<-c(NA,NA);dif<-c(NA,NA)
for (i in 1:nrow(data)){
  p[i]<-wilcox.test(as.numeric(afr[i,]),as.numeric(eur[i,]),na.rm=T)$p.value
  dif[i]<-mean(as.numeric(eur[i,]),na.rm=T)-mean(as.numeric(afr[i,]),na.rm=T)
}
fdr<-p.adjust(p,method="fdr")
re<-data.frame(rownames(data),dif,p,fdr)
write.table(re,file="re_unadjusted.txt",sep="\t",quote=F)



re<-read.table(gzfile("re_unadjusted.txt.gz"))
target<-read.table("../stableCpG_inWGBS.txt")
re1<-re[re$rownames.data. %in% target$rank,]
summary(abs(re1$dif))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0004602 0.0079368 0.0180324 0.0217789 0.0320545 0.0830463
re$fdr
# [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
# [47] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
summary(re1$p)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.02462 0.26645 0.39770 0.47127 0.71508 1.00000 
summary(abs(re$dif))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.01963 0.04357 0.05963 0.08314 0.92691
summary(abs(re$dif[re$fdr<0.05]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0264  0.3031  0.5108  0.4702  0.6513  0.9269    1281 

all<-re
# all$dif<-abs(all$dif)
all$group[all$rownames.data. %in% target$rank]<-"Non-pop-specific"
# all$group[all$fdr<0.05]<-"Pop-specific"
all$group[is.na(all$group)]<-"Other"

g<-ggplot(all, aes(x=dif, color=group)) +
  geom_density()+theme_classic()+scale_color_manual(values=c("#4DBBD5","#999999"))+xlab("DNAm changes between the two populations")
ggsave(g,file="density_new.pdf",width = 5,height = 3)






