meth<-read.table(gzfile("../../Old_CMR_analysis/data_CMR.txt.gz"))
tran<-read.table(gzfile("../../../../MHB_May3/WGBS_transcript/combat_total_WGBSandArray.txt.gz"),header = T,row.names = 1)
colnames(tran)<-sub(".cel.gz", "", colnames(tran))
colnames(tran)<-sub(".CEL.gz", "", colnames(tran))
colnames(tran)<-sub("_b.*_lcl_hg133plus2", "", colnames(tran))
pheno<-read.table("README.txt",sep="\t")
meth1<-meth[,pheno$V1]
tran1<-tran[,pheno$V3]
score<-read.table(gzfile("../../Old_CMR_analysis/PCAscore.txt.gz"))
score1<-score[pheno$V1,]
pair<-read.table("gene-meth_all.txt")
tran[tran$symbol %in% pair$V2,1:2]
pair<-read.table("pair_all.txt")

p<-c(NA,NA);dif<-c(NA,NA);r<-c(NA,NA)
for (i in 1:nrow(pair)){
  res<-lm(as.numeric(meth1[pair$V4[i],-c(2)])~as.numeric(tran1[pair$V1[i],-c(2)])+as.numeric(pheno$V2[-c(2)])+pheno$V4[-c(2)]+score1[-c(2),1])#+score1[,2]+score1[,3])
  # res<-lm(as.numeric(tran1[pair$V1[i],])~as.numeric(meth1[pair$V4[i],])+pheno$V2+pheno$V4)#+score1[,1]+score1[,2]+score1[,3])
  dif[i]<-summary(res)$coefficients[2,1]
  p[i]<-summary(res)$coefficients[2,4]
  r[i]<-cor(as.numeric(tran1[pair$V1[i],]),as.numeric(meth1[pair$V4[i],]))
  # dif[i]<-cor.test(as.numeric(tran1[pair$V1[i],-c(2)]),as.numeric(meth1[pair$V4[i],-c(2)]),method="spearman")$estimate
  # p[i]<-cor.test(as.numeric(tran1[pair$V1[i],]),as.numeric(meth1[pair$V4[i],]),method="spearman")$p.value
}
fdr<-p.adjust(p,method="fdr")
re<-data.frame(pair,dif,p,fdr,r)
re[re$p<0.05 & abs(re$r)>0.3,]
# V1              V2    V3                    V4         dif          p       fdr          r
# 20 211515_s_at ENSG00000133059 DSTYK 1:205141231:205141245  1.03256625 0.01657375 0.6602026  0.6020571
# 35 211067_s_at ENSG00000007237  GAS7  17:10056902:10057753  0.17444047 0.03441987 0.6602026  0.7549465
# 61    41644_at ENSG00000111961 SASH1 6:148663571:148663585 -0.08030816 0.03242352 0.6602026 -0.6480752
# 72   212770_at ENSG00000140332  TLE3  15:70387017:70387094 -0.43147322 0.03927174 0.6602026 -0.3267740

# p<-c(NA,NA);dif<-c(NA,NA)
# for (i in 1:nrow(pair)){
#   res<-lm(as.numeric(meth1[pair$V4[i],1:5])~as.numeric(tran1[pair$V1[i],1:5])+pheno$V2[1:5]+pheno$V4[1:5])#+score1[,1]+score1[,2]+score1[,3])
#   # res<-lm(as.numeric(tran1[pair$V1[i],1:5])~as.numeric(meth1[pair$V4[i],1:5]))#+pheno$V2[1:5])
#   dif[i]<-summary(res)$coefficients[2,1]
#   p[i]<-summary(res)$coefficients[2,4]
#   # p[i]<-cor.test(as.numeric(tran1[pair$V1[i],1:5]),as.numeric(meth1[pair$V4[i],1:5]))$p.value
# }
# fdr<-p.adjust(p,method="fdr")
# re<-data.frame(pair,dif,p,fdr)
# re[re$fdr<0.05,]
# # None
# 
# p<-c(NA,NA);dif<-c(NA,NA)
# for (i in 1:nrow(pair)){
#   # res<-lm(as.numeric(meth1[pair$V4[i],])~as.numeric(tran1[pair$V1[i],])+pheno$V2+pheno$V4)#+score1[,1]+score1[,2]+score1[,3])
#   res<-lm(as.numeric(tran1[pair$V1[i],6:11])~as.numeric(meth1[pair$V4[i],6:11])+pheno$V2[6:11])
#   dif[i]<-summary(res)$coefficients[2,1]
#   p[i]<-summary(res)$coefficients[2,4]
#   # p[i]<-cor.test(as.numeric(tran1[pair$V1[i],6:11]),as.numeric(meth1[pair$V4[i],6:11]))$p.value
# }
# fdr<-p.adjust(p,method="fdr")
# re<-data.frame(pair,dif,p,fdr)
# re[re$fdr<0.05,]
# None

info_expression<-read.table(gzfile("../../../../MHB_May3/WGBS_transcript/pheno.txt.gz"))
rownames(info_expression)<-sub(".cel.gz", "", rownames(info_expression))
rownames(info_expression)<-sub(".CEL.gz", "", rownames(info_expression))
rownames(info_expression)<-sub("_b.*_lcl_hg133plus2", "", rownames(info_expression))
expression_AFR<-tran[,rownames(info_expression[info_expression$subgroups==1,])] # 156
expression_EUR<-tran[,rownames(info_expression[info_expression$subgroups==0,])] # 149

expression_AFR<-expression_AFR[unique(pair$V1),]
expression_EUR<-expression_EUR[unique(pair$V1),] # -which(names(expression_EUR) %in% "GSM1202824")
p<-c(0,0);
for (i in 1:(nrow(expression_AFR))) {
  #p[i]<-wilcox.test(jitter(as.numeric(expression_AFR[i,])),jitter(as.numeric(expression_EUR[i,])))$p.value
  p[i]<-wilcox.test(as.numeric(expression_AFR[i,]),as.numeric(expression_EUR[i,]))$p.value
  
}
fdr<-p.adjust(p,method = "fdr")
re_w<-data.frame(unique(pair$V1),p,fdr)

write.table(re_w,file = "result_expressiondif_all.txt",quote = F,row.names = F,sep="\t")

expression_AFR<-expression_AFR[re$V1[re$p<0.05 & abs(re$r)>0.3],]
expression_EUR<-expression_EUR[re$V1[re$p<0.05 & abs(re$r)>0.3],] # -which(names(expression_EUR) %in% "GSM1202824")
p<-c(0,0);
for (i in 1:(nrow(expression_AFR))) {
  #p[i]<-wilcox.test(jitter(as.numeric(expression_AFR[i,])),jitter(as.numeric(expression_EUR[i,])))$p.value
  p[i]<-wilcox.test(as.numeric(expression_AFR[i,]),as.numeric(expression_EUR[i,]))$p.value
  
}
fdr<-p.adjust(p,method = "fdr")
re_w<-data.frame(re$V1[re$p<0.05 & abs(re$r)>0.3],p,fdr)

write.table(re_w,file = "result_expressiondif_sig.txt",quote = F,row.names = F,sep="\t")

re_w<-re_w[re_w$fdr<0.05,]
# re.V1.re.p...0.05...abs.re.r....0.3.           p         fdr
# 2                          211067_s_at 0.002533725 0.005067451
# 3                             41644_at 0.002512083 0.005067451


library(gplots)
library(ggpubr)
probenames<-re_w$re.V1.re.p...0.05.
gene<-c("PSPH","SCARF1","PIEZO1")
for (i in 1:(length(probenames))) {
  afr <- data.frame(as.numeric(expression_AFR[probenames[i],]),"AFR")
  names(afr)<-c("Expression","Race")
  eur <- data.frame(as.numeric(expression_EUR[probenames[i],]),"EUR")
  names(eur)<-c("Expression","Race")
  df<-rbind(afr,eur)
  p <- ggboxplot(df, x = "Race", y = "Expression",
                 color = "Race", palette =c('#E69F00', '#009E73'),show.legend=FALSE,
                 add = "jitter")+ylab(paste("Expression levels for",gene[i]))+
    theme(axis.title.x=element_blank(),
          legend.position = "none")#+ stat_compare_means()
  ggsave(plot = p, width = 4, height = 4, dpi = 300, file=paste(probenames[i], ".pdf", sep=""))
}

