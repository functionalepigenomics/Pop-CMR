data<-read.table(gzfile("data_CMR.txt.gz"))
# data<-data[,order(names(data))]
pheno<-read.table(gzfile("README.txt.gz"),header = T)
rownames(pheno)<-paste(pheno$race,pheno$Sample,sep=".")
pheno<-pheno[names(data),]
# pheno$sex[pheno$sex=="F"]<-0
# pheno$sex[pheno$sex=="M"]<-1
pheno$race[pheno$race=="black"]<-0
pheno$race[pheno$race=="white"]<-1
pheno$sex<-as.numeric(pheno$sex)
pheno$race<-as.numeric(pheno$race)

score<-read.table("PCAscore.txt")
score<-score[names(data),]
cor<-c(NA,NA);p<-c(NA,NA)
for (i in 1:13){
  re<-cor.test(score[,i],pheno[,2])
  p[i]<-re$p.value
  cor[i]<-re$estimate
}
# > p
# [1] 0.69867924 0.15358365 0.73103046 0.18105812 0.39909905 0.01631986 0.17706476
# [8] 0.39983996 0.44041660 0.92331154 0.55025619 0.98929934 0.36393038
# > cor
# [1] -0.08532771 -0.30742886 -0.07580186  0.28901076  0.18460202 -0.49504218 -0.29156139
# [8] -0.18431711  0.16913460 -0.02125542 -0.13134460  0.00296170  0.19848611

cor<-c(NA,NA);p<-c(NA,NA)
for (i in 1:13){
  re<-cor.test(score[,i],pheno[,3])
  p[i]<-re$p.value
  cor[i]<-re$estimate
}
# > p
# [1] 0.001898393 0.962001254 0.578974199 0.893258768 0.456953782 0.996275765 0.770608115
# [8] 0.853255774 0.920952572 0.994819447 0.122079873 0.484192823 0.371935553
# > cor
# [1] -0.52126950  0.05913911  0.12510089  0.08439811 -0.14841793 -0.04710563  0.11156783
# [8] -0.12581840 -0.01419993  0.12280603 -0.26965907  0.30139168 -0.30048914

cor<-c(NA,NA);p<-c(NA,NA)
for (i in 1:13){
  re<-cor.test(score[,i],pheno[,4])
  p[i]<-re$p.value
  cor[i]<-re$estimate
}
# > p
# [1] 0.31692096 0.63697440 0.56742709 0.21343828 0.34607911 0.23787233 0.52719907
# [8] 0.10899679 0.80569923 0.50849205 0.03827761 0.29565134 0.87990418
# > cor
# [1]  0.21832603 -0.10393367 -0.12577685 -0.26962533 -0.20583037 -0.25626861 -0.13894768
# [8] -0.34308797 -0.05427919  0.14523235 -0.43451646  0.22788779 -0.03335607

cor<-c(NA,NA);p<-c(NA,NA)
for (i in 1:13){
  re<-cor.test(score[,i],pheno[,5])
  p[i]<-re$p.value
  cor[i]<-re$estimate
}
# > p
# [1] 0.48944465 0.21433959 0.63032825 0.77851723 0.53005574 0.13795984 0.01620889
# [8] 0.89352408 0.60673218 0.45473079 0.81549941 0.34143931 0.27710719
# > cor
# [1]  0.15174792 -0.26911509 -0.10597724  0.06204883  0.13799739 -0.31895843 -0.49548726
# [8] -0.02955084  0.11330329  0.16395915  0.05149364 -0.20777552 -0.23657870

cor<-c(NA,NA);p<-c(NA,NA)
for (i in 1:13){
  re<-cor.test(score[,i],pheno[,6])
  p[i]<-re$p.value
  cor[i]<-re$estimate
}
# > p
# [1] 9.245797e-12 9.865223e-01 4.853734e-01 7.641743e-01 9.956661e-01 4.703360e-01
# [7] 8.092551e-01 7.803370e-01 5.407378e-01 8.951954e-01 9.199092e-01 9.764670e-01
# [13] 6.839743e-01
# > cor
# [1]  0.946214892 -0.003730370 -0.153156778  0.066176994  0.001199489  0.158412867
# [7] -0.053267572 -0.061526528  0.134465023 -0.029084482 -0.022201312 -0.006514106
# [13] -0.089706034

dif<-c(NA,NA);p<-c(NA,NA);se<-c(NA,NA);
for (i in 1:nrow(data)){
  res<-lm(as.numeric(data[i,])~pheno$race +
            pheno$sex+pheno$age)#+score[,2]+score[,3])#+pheno$height+pheno$weight)#+
  #   +score[,5]+score[,6]+score[,7]+score[,8]+
  # score[,9]+score[,10]+score[,11]+score[,12])
  dif[i]<-summary(res)$coefficients[2,1]
  # dif[i]<-mean(as.numeric(data[i,7:23]))-mean(as.numeric(data[i,1:6]))
  # p[i]<-wilcox.test(as.numeric(data[i,1:6]),as.numeric(data[i,7:23]))$p.value
  p[i]<-summary(res)$coefficients[2,4]
  se[i]<-summary(res)$coefficients[2,2]
}
fdr<-p.adjust(p,method="fdr")
re<-data.frame(dif,se,p,fdr)
rownames(re)<-rownames(data)
nrow(re[re$fdr<0.05,]) # 26
nrow(re[re$fdr<0.05 & abs(re$dif) > 0.05,]) # 24
nrow(re[re$fdr<0.05 & abs(re$dif) > 0.1,]) # 21
# write.table(re,file="re_WithoutCMRvalidate.txt",sep="\t",quote = F)
write.table(re,file="re_onlysexandage.txt",sep="\t",quote = F)


cmr<-read.table("CMR_cfDNA.bed")
data<-data[unique(cmr$V1),]
dif<-c(NA,NA);p<-c(NA,NA);se<-c(NA,NA);
for (i in 1:nrow(data)){
  res<-lm(as.numeric(data[i,])~pheno$race +
            pheno$sex+pheno$age+score[,2]+score[,3])#+pheno$height+pheno$weight)#+
          #   +score[,5]+score[,6]+score[,7]+score[,8]+
          # score[,9]+score[,10]+score[,11]+score[,12])
  dif[i]<-summary(res)$coefficients[2,1]
  # dif[i]<-mean(as.numeric(data[i,7:23]))-mean(as.numeric(data[i,1:6]))
  # p[i]<-wilcox.test(as.numeric(data[i,1:6]),as.numeric(data[i,7:23]))$p.value
  p[i]<-summary(res)$coefficients[2,4]
  se[i]<-summary(res)$coefficients[2,2]
}
fdr<-p.adjust(p,method="fdr")
re<-data.frame(dif,se,p,fdr)
rownames(re)<-rownames(data)
nrow(re[re$fdr<0.05,]) # 22
nrow(re[re$fdr<0.05 & abs(re$dif) > 0.05,]) # 22
nrow(re[re$fdr<0.05 & abs(re$dif) > 0.1,]) # 20
write.table(re,file="re_CMRvalidate.txt",sep="\t",quote = F)