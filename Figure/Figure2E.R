setwd("~/KoborLab/kobor_space/zdong/Population_dataset/MHB_New/Analysis/CMR_analysis/Permutation_PCA")
library(gplots) 
library(devtools)
library(ggbiplot)
test<-read.table(gzfile("../../Old_CMR_analysis/data_CMR.txt.gz"),row.names=1,header=T)
test<-test[,1:62] ##only for MHBmap data

test1<-test;
a<-c(0,0);b<-c(0,0)
####========================== PCA =============================
for (i in 1:1000){
  test<-test1[sample(nrow(test1), 101), ]
case<-test[,c(55:62)]
control<-test[,c(1:54)]
case_pca=rbind(case,c(1,1))
control_pca=rbind(control,c(0,0))
test_pca<-data.frame(case_pca,control_pca)
#test_pca<-as.matrix(t(test_pca))
rownames(test_pca)[length(rownames(test_pca))]<-"Population"
iris<-as.data.frame(t(test_pca))

#!!!!Imputation for PCA
ir <- (iris[, 1:(ncol(iris)-1)])
if (sum(is.na(ir))>0){
  cM <- colMeans(ir, na.rm=TRUE)
  indx <- which(is.na(ir), arr.ind=TRUE)
  ir[indx] <- cM[indx[,2]]}
ir[mapply(is.infinite, ir)] <--99999
# log transform 
log.ir=log(ir)
#log.ir[mapply(is.infinite, log.ir)] <--99999
ir.species <- iris[, ncol(iris)]
ir.species[ir.species==1]<-"EUR"
ir.species[ir.species==0]<-"AFR"
# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 
ir.pca <- prcomp(ir,
                 center = TRUE,
                 scale. = TRUE)
#summary(ir.pca)$importance[,1:6]
a[i]<-summary(ir.pca)$importance[2,1]
b[i]<-summary(ir.pca)$importance[2,2]
}

n<-1:1000
c<-a+b
final<-data.frame(n,a,b,c)


summary(a)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.08871 0.11730 0.12445 0.12519 0.13285 0.16398 
wilcox.test(a,mu=0.315,alternative = "less")$p.value # 1.665307e-165

theme_set(theme_classic())
g <- ggplot(final,aes(a))+
  geom_density(fill='#636363',color="#636363",alpha = 0.8#,position = "stack"
  )+ labs(x="Proportion of variance explained by the top one PC")+
  # geom_rug(aes(x = cor, y = 0,colour=factor(num)),
  #          position = position_jitter(height = 0))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y  = element_text(size=12),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size=11),
        legend.position=c(0.12,0.9),legend.title =  element_blank(),
        legend.spacing.x = unit(0.1, 'cm'),
        axis.text.y = element_text(size=12),
        axis.text.x =element_text(size=12),
        axis.ticks.x = element_line(colour = "black", size = 0.8),
        axis.ticks.y=element_line(colour = "black", size = 0.8),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())+scale_fill_manual(values=c('#636363'))#+
# geom_vline(aes(xintercept=0.3),color="black", linetype="dashed", size=0.8)
# scale_color_manual(values=c("#E69F00","#009E73"))
ggsave(plot = g, width = 4,height = 3.4, dpi = 300, filename = "density.pdf")


