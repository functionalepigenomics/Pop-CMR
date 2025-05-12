setwd("~/KoborLab/kobor_space/zdong/Population_dataset/MHB_New/Analysis/CMR_analysis")
library(gplots) 
test<-read.table(gzfile("../Old_CMR_analysis/data_CMR.txt.gz"))
test<-test[,1:62]
sig<-read.table("DM-CMR.txt.gz")
test<-test[rownames(sig),]

# http://www.rapidtables.com/web/color/RGB_Color.htm
condition_colors <- unlist(lapply(colnames(test),function(x){
  if(grepl('AFR',x)) '#E69F00' #pink
  else if(grepl('EUR',x)) '#009E73' #grey
  #else if(grepl('ceu',x)) '#008000'
  #else if(grepl('H.s',x)) '#0000FF'
  # else '#FF0000'
}))
input<-as.matrix(test)
myCol <- colorRampPalette(c("blue","yellow"))(100)
myBreaks <- c(0,0.2,0.4,0.6,0.8,1)
pdf("test_aml0.1.pdf")
heatmap.2(input,na.rm=T,#Colv=NA,Rowv=NA,
          #colCol=c("#CC0066","#00994C"),
          labRow=NA,lhei = c(1.09,5),lwid = c(1,4),key.title=NA,key.ylab=NA,
          labCol=NA,
          trace="none", density="none", col=myCol,symbreaks = F,cexRow=1, cexCol=0.2, margins = c(1,1),
          #Euclidean distance with Ward's linkage
          distfun=function(x) dist(x, method="euclidean"), hclustfun=function(x) hclust(x, method="ward.D2"),#ward.D2
          #1 minus Pearson correlation distance with average linkage
          #distfun=function(x) as.dist(1-cor(t(x))), hclustfun=function(x) hclust(x, method="average"),
          keysize=1.1,key.xlab="AML",ColSideColors=condition_colors, scale="none")#,hclust=function(x) hclust(x,method="average"),distfun=function(x) as.dist((1-cor(x))/2))
legend(0.09,0.9421,bty = "n",text.width=c(0.00,0.429),legend=c("EUR","AFR"),inset=c(-0.05,-0.1),lty=0,xpd=T,horiz=T,cex=0.85)
dev.off()

####========================== PCA =============================
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

# plot method
#plot(ir.pca, type = "l")
# summary method
summary(ir.pca)$importance[,1:6]

regions<-c("Utah","Utah","Stockholm","Stockholm","Stockholm","Stockholm","Ibadan","Ibadan",
           "Ibadan",'Ibadan','Ibadan',"Ibadan",'Ibadan','Ibadan',"Ibadan",'Ibadan')
library(devtools)
library(ggbiplot)
#pdf("pca.pdf")
theme_set(theme_classic())
g<-ggbiplot(ir.pca, obs.scale = 1, var.scale = 1,#labels = regions,
            groups = ir.species, ellipse = TRUE,  var.axes = F,
            circle = TRUE) +
  #scale_color_discrete(name ='') +
  theme(legend.direction = 'horizontal',
        legend.position = 'top')+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),#element_rect(colour = "black", size=2),
        # axis.title = element_text(size=12),legend.text = element_text(size=12),
        legend.position = c(0.48,1.12),axis.ticks = element_line(size = 1.2),
        legend.title = element_text(size=12),
        axis.text =element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12)
        )+
  scale_colour_manual(values=c("#E69F00","#009E73"))+
  labs(color="Population",size=12)+ coord_fixed(ratio=1)
#dev.off()
ggsave(plot = g, width = 2.7, height = 2.7, dpi = 300, filename = "pca_aml0.1.pdf")
#ggsave(plot = g, width = 5, height = 5, dpi = 300, filename = "pca_new_avo_regions.pdf")

# ####========================== PCA _ region =============================
# case<-test[,c(11:16)]
# control<-test[,c(1:10)]
# case_pca=rbind(case,c(1,1))
# control_pca=rbind(control,c(0,0))
# test_pca<-data.frame(case_pca,control_pca)
# #test_pca<-as.matrix(t(test_pca))
# rownames(test_pca)[length(rownames(test_pca))]<-"Population"
# iris<-as.data.frame(t(test_pca))
# 
# #!!!!Imputation for PCA
# ir <- (iris[, 1:(ncol(iris)-1)])
# if (sum(is.na(ir))>0){
#   cM <- colMeans(ir, na.rm=TRUE)
#   indx <- which(is.na(ir), arr.ind=TRUE)
#   ir[indx] <- cM[indx[,2]]}
# ir[mapply(is.infinite, ir)] <--99999
# # log transform 
# log.ir=log(ir)
# #log.ir[mapply(is.infinite, log.ir)] <--99999
# ir.species <- iris[, ncol(iris)]
# ir.species[ir.species==1]<-"EUR"
# ir.species[ir.species==0]<-"AFR"
# # apply PCA - scale. = TRUE is highly 
# # advisable, but default is FALSE. 
# ir.pca <- prcomp(ir,
#                  center = TRUE,
#                  scale. = TRUE)
# 
# # plot method
# #plot(ir.pca, type = "l")
# # summary method
# summary(ir.pca)$importance[,1:6]
# 
# regions<-c("Utah","Utah","Stockholm","Stockholm","Stockholm","Stockholm","Ibadan","Ibadan",
#            "Ibadan",'Ibadan','Ibadan',"Ibadan",'Ibadan','Ibadan',"Ibadan",'Ibadan')
# library(devtools)
# library(ggbiplot)
# #pdf("pca.pdf")
# theme_set(theme_classic())
# g<-ggbiplot(ir.pca, obs.scale = 1, var.scale = 1,labels = regions,
#             groups = ir.species, ellipse = TRUE,  var.axes = F,
#             circle = TRUE) +
#   #scale_color_discrete(name ='') +
#   theme(legend.direction = 'horizontal',
#         legend.position = 'top')+
#   theme(panel.border = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_rect(colour = "black", size=2),
#         axis.title = element_text(size=12),legend.text = element_text(size=12),
#         legend.position = c(0.48,1.12),axis.ticks = element_line(size = 1.2),
#         #legend.title = element_text(size=12),
#         axis.text =element_text(size=12),
#         axis.line.x = element_blank(),
#         axis.line.y = element_blank())+
#   scale_colour_manual(values=c("#E69F00","#009E73"))+
#   labs(color="Population",size=12)
# #dev.off()
# #ggsave(plot = g, width = 5, height = 5, dpi = 300, filename = "pca_new_avo_2.pdf")
# ggsave(plot = g, width = 5, height = 5, dpi = 300, filename = "pca_new_avo_regions_beta0.1.pdf")
