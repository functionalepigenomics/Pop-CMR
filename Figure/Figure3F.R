geno<-read.table("ExtractSNP_52.raw",header=T)
all<-read.table("../../../Old_CMR_analysis/Genotype/Allale/3202.txt",header=T,sep="\t")
unrelated<-read.table("../../../Old_CMR_analysis/Genotype/Allale/1000G_2504_high_coverage.sequence.index",header=T)
sel<-all[all$Sample.name %in% unrelated$SAMPLE_NAME,]

eur<-geno[geno$IID %in% sel$Sample.name[sel$Superpopulation.code=="EUR"],]
afr<-geno[geno$IID %in% sel$Sample.name[sel$Superpopulation.code=="AFR"],]
write.table(eur$IID,file="sample_eur.txt",row.names = F,col.names = F,quote=F)
write.table(afr$IID,file="sample_afr.txt",row.names = F,col.names = F,quote=F)


p<-c(NA,NA);allele.fre.afr<-c(NA,NA);allele.fre.eur<-c(NA,NA);
for (i in 7:ncol(geno)){
  allele.fre.afr[i-6]<-sum(afr[,i])/(nrow(afr)*2)
  allele.fre.eur[i-6]<-sum(eur[,i])/(nrow(eur)*2)
  dd<-c(sum(eur[,i]),(nrow(eur)*2),sum(afr[,i]),(nrow(afr)*2))
  dim(dd)<-c(2,2)
  p[i-6]<-fisher.test(dd)$p.value
}
fdr<-p.adjust(p,method="fdr")
re<-data.frame(names(geno)[7:ncol(geno)],allele.fre.eur,allele.fre.afr,p,fdr)
re$dif<-abs(re$allele.fre.eur-re$allele.fre.afr)
nrow(re[re$fdr<0.05 & re$dif > 0.1,])/nrow(re) # 0.6923077 # 36
write.table(re,file="re.txt",sep="\t",quote = F,row.names = F)

geno1<-read.table("ExtractSNP.raw",header=T)
eur1<-geno1[geno1$IID %in% sel$Sample.name[sel$Superpopulation.code=="EUR"],]
afr1<-geno1[geno1$IID %in% sel$Sample.name[sel$Superpopulation.code=="AFR"],]
num<-c(NA,NA);n=0
for (j in seq(7,ncol(geno1),by=52)){n=n+1;
afr<-afr1[,j:(j+51)]
eur<-eur1[,j:(j+51)]
p<-c(NA,NA);allele.fre.afr<-c(NA,NA);allele.fre.eur<-c(NA,NA);
for (i in 1:ncol(afr)){
  allele.fre.afr[i]<-sum(afr[,i])/(nrow(afr)*2)
  allele.fre.eur[i]<-sum(eur[,i])/(nrow(eur)*2)
  dd<-c(sum(eur[,i]),(nrow(eur)*2),sum(afr[,i]),(nrow(afr)*2))
  dim(dd)<-c(2,2)
  p[i]<-fisher.test(dd)$p.value
}
fdr<-p.adjust(p,method="fdr")
re<-data.frame(allele.fre.eur,allele.fre.afr,p,fdr)
re$dif<-abs(re$allele.fre.eur-re$allele.fre.afr)
num[n]<-length(re$fdr[re$fdr<0.05 & re$dif > 0.1])
# num[n]<-length(re$fdr[re$fdr<5e-8])
}
mean(num) # 18.2
summary(num) # re$fdr<0.05 & re$dif > 0.1
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.0    16.0    18.0    18.2    21.0    29.0 
mean(num)/52 # 0.35
length(num[num>=36])/length(num) # 0

#================================== bar chart ============================================
library(ggplot2)
dfm<-data.frame(name=c("SNPs for Pop-CMRs","Random SNPs"),mpg=c(0.6923077,0.35),cyl=c("SNPs for Pop-CMRs","Random SNPs"))
dfm$name <- factor(dfm$name, levels = c("SNPs for Pop-CMRs","Random SNPs"))
p<-ggplot(data=dfm, aes(x=name,y=mpg))  +theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        axis.text = element_text(size=12),
        legend.title =element_blank()
  )+ylab("Proportion of population-specific loci")+xlab("")+
  #geom_rect(aes(),xmin =-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,alpha = 0.1) +
  geom_bar(aes(fill = name),stat = "identity",position = "dodge") + ylim(0,1)+
  scale_fill_manual(values=c("#e41a1c","#fc9272"))#+
#facet_wrap(~ tissue, nrow=2)
ggsave("WGBS.pop-allele.bbar.chart.pdf",p, width=4, height=2.6, dpi = 300)