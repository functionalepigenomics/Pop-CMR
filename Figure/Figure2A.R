score<-read.table(gzfile("../Old_CMR_analysis/PCAscore.txt.gz"))
# PC1      PC2      PC3      PC4      PC5      PC6
# Standard deviation     65.32411 42.07543 36.18749 32.68246 31.81712 30.33773
# Proportion of Variance  0.11641  0.04829  0.03572  0.02914  0.02762  0.02511
# Cumulative Proportion   0.11641  0.16470  0.20043  0.22957  0.25718  0.28229
pheno<-read.table(gzfile("../Old_CMR_analysis/README.txt.gz"),sep="\t",header = T)
rownames(pheno)<-pheno$Sample
pheno<-pheno[rownames(score),]
cor<-c(NA,NA);p<-c(NA,NA)
for (i in 1:10){
  re<-cor.test(score[,i],pheno[,2])
  p[i]<-re$p.value
  cor[i]<-re$estimate
}

# > p
# [1] 0.4760934 0.3332849 0.1771931 0.4954477 0.3668654 0.4108412 0.2519910 0.6239908 0.8006441
# [10] 0.8184508
# > cor
# [1] -0.09218333 -0.12493522  0.17360854  0.08819814 -0.11658101  0.10631381  0.14768817
# [8] -0.06348486  0.03272669 -0.02975035
cor1<-c(NA,NA);p1<-c(NA,NA)
for (i in 1:10){
  re<-cor.test(score[,i],pheno[,3])
  p1[i]<-re$p.value
  cor1[i]<-re$estimate
}
# > cor1
# [1] -0.18629878  0.20721457 -0.28557141 -0.02775735  0.20240789  0.26096505 -0.07741967
# [8]  0.15827816  0.07818212  0.10528314
# > p1
# [1] 0.14711673 0.10609488 0.02445485 0.83042628 0.11463327 0.04049430 0.54977599 0.21919073
# [9] 0.54584093 0.41541801

data<-read.table(gzfile("../Old_CMR_analysis/data_CMR.txt.gz"))
data<-data[,rownames(score)]
dif<-c(NA,NA);p<-c(NA,NA);se<-c(NA,NA);
for (i in 1:nrow(data)){
  res<-lm(as.numeric(data[i,])~pheno[,3]+pheno[,2]+score[,1])
  dif[i]<-summary(res)$coefficients[2,1]
  p[i]<-summary(res)$coefficients[2,4]
  se[i]<-summary(res)$coefficients[2,2]
}
fdr<-p.adjust(p,method="fdr")
re<-data.frame(dif,se,p,fdr)
rownames(re)<-rownames(data)
write.table(re,file="re_PC1.txt",sep="\t",quote = F)


nrow(re[abs(re$dif)>0.03 & re$fdr<0.05,]) # 137
nrow(re[abs(re$dif)>0.05 & re$fdr<0.05,]) # 136
nrow(re[abs(re$dif)>0.1 & re$fdr<0.05,]) # 135
sig<-re[abs(re$dif)>0.1 & re$fdr<0.05,]
write.table(sig,file="re_aml0.1andfdr0.05_PC1.txt",sep="\t",quote = F)

### permutation
a2<-"na";a2p<-"na"
new_data<-data[rownames(sig),]
#data1<-replace(new_data,is.na(new_data),0)
afr=new_data[,9:62];
eur=new_data[,1:8];
for (j in 1:1000){a<-c(1,2);
d<-afr[,sample(ncol(afr),ncol(eur))];
d<-data.frame(eur,d)
p<-c(NA,NA);dif<-c(NA,NA)
for (i in 1:nrow(afr)){
  res<-lm(as.numeric(d[i,])~pheno[names(d),3]+pheno[names(d),2]+score[names(d),1])
  dif[i]<-summary(res)$coefficients[2,1]
  p[i]<-summary(res)$coefficients[2,4]
  # se[i]<-summary(res)$coefficients[2,2]
}
a2<-data.frame(a2,dif)
fdr<-p.adjust(p,method="fdr")
a2p<-data.frame(a2p,fdr)
}
a3<-a2[,2:1001]
a3p<-a2p[,2:1001]
require(matrixStats)
dif_median<-matrixStats::rowMedians(as.matrix(a3))
length(dif_median[dif_median > 0.1 | dif_median < -0.1])/length(dif_median) # 1
dif_low<-matrixStats::rowMins(abs(as.matrix(a3)))
length(dif_low[dif_low>0.1])/length(dif_low) # 0.5333333
length(dif_low[dif_low>0.05])/length(dif_low) # 0.7925926
length(dif_low[dif_low>0.03])/length(dif_low) # 0.8666667

fdr_median<-matrixStats::rowMedians(as.matrix(a3p))
length(fdr_median[fdr_median<0.05])/length(fdr_median) # 0.7555556
fdr_max<-matrixStats::rowMaxs(abs(as.matrix(a3p)))
length(fdr_max[fdr_max<0.05])/length(fdr_max) # 0.04444444

per<-data.frame(rownames(sig),dif_median,fdr_median)
re_per<-per[abs(per$dif_median)>0.1 & per$fdr_median<0.05,]

write.table(re_per,file="re_aml0.1andfdr0.05_PC1_permutation.txt",sep="\t",quote = F)

eur=new_data[re_per$rownames.sig.,1:8];
sub<-c(1,1,1,1,1,0,0,0)
p<-c(NA,NA);dif<-c(NA,NA)
for (i in 1:nrow(eur)){
  res<-lm(as.numeric(eur[i,])~sub+pheno[names(eur),2]+score[names(eur),1])
  dif[i]<-summary(res)$coefficients[2,1]
  p[i]<-summary(res)$coefficients[2,4]
  # se[i]<-summary(res)$coefficients[2,2]
}
eur<-p.adjust(p,method="fdr")
eur<-data.frame(re_per$rownames.sig.,dif,fdr)
re_eur<-eur[abs(afr$dif)>0.1 & eur$fdr<0.05,] # 0

## permutation in afr
afr=new_data[re_per$rownames.sig.,9:62]
gam<-afr[,8:53]
yri<-afr[,c(1:7,54)]
a2<-"na";a2p<-"na"
for (j in 1:1000){a<-c(1,2);
d<-gam[,sample(ncol(gam),ncol(yri))];
d<-data.frame(yri,d)
p<-c(NA,NA);dif<-c(NA,NA)
for (i in 1:nrow(yri)){
  res<-lm(as.numeric(d[i,])~c(replicate(ncol(yri),1),replicate(ncol(yri),0))+pheno[names(d),2]+score[names(d),1])
  dif[i]<-summary(res)$coefficients[2,1]
  p[i]<-summary(res)$coefficients[2,4]
  # se[i]<-summary(res)$coefficients[2,2]
}
a2<-data.frame(a2,dif)
fdr<-p.adjust(p,method="fdr")
a2p<-data.frame(a2p,fdr)
}
a3<-a2[,2:1001]
a3p<-a2p[,2:1001]
require(matrixStats)
dif_median<-matrixStats::rowMedians(as.matrix(a3))
# dif_low<-matrixStats::rowMins(abs(as.matrix(a3)))
fdr_median<-matrixStats::rowMedians(as.matrix(a3p))
per<-data.frame(re_per$rownames.sig.,dif_median,fdr_median)
re_per<-per[abs(per$dif_median)>0.1 & per$fdr_median<0.05,]
re_per
# re_per.rownames.sig. dif_median  fdr_median
# 81 10:12171054:12171110 -0.2333084 0.002604955

## check with previous 64 sites
re<-read.table("re_aml0.1andfdr0.05_PC1.txt.gz")
re_per<-read.table("re_aml0.1andfdr0.05_PC1_permutation.txt.gz")
re<-re[re_per$rownames.sig.,]
final<-re[-which(rownames(re) %in% c("10:12171054:12171110")),]
write.table(final,file="DM-CMR.txt",sep="\t",quote = F)


pre<-read.table("../Old_CMR_analysis/re_aml0.1andfdr0.05.txt.gz")
length(intersect(rownames(pre),rownames(final))) # 43

### pca analysis
library(ggrepel)
re<-read.table("re_PC1.txt")
# add a column of NAs
re$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
re$diffexpressed[re$dif > 0.1 & re$fdr < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
re$diffexpressed[re$dif < -0.1 & re$fdr < 0.05] <- "DOWN"
re$diffexpressed[-which(rownames(re) %in% rownames(final))] <- "NO"
# re$delabel <- NA
# re$delabel[re$diffexpressed != "NO"] <- re$cpg[re$diffexpressed != "NO"]
g<-ggplot(data=re, aes(x=dif, y=-log10(p), col=diffexpressed)) +
  geom_point() + 
  theme_minimal() +
  #geom_text_repel() +
  scale_color_manual(values=c("#E64B35", "#bdbdbd", "red")) +
  geom_vline(xintercept=c(-0.1, 0.1), col="black",linetype ="dashed") +
  geom_hline(yintercept=-log10(0.000182100402994445), col="black",linetype ="dashed")+
  theme(legend.position = "none",
        axis.text = element_text(size = 11)
  )+xlab("AML change (EUR-AFR)")
ggsave(g,filename = "valcano.pdf",width = 4,height = 3.4)