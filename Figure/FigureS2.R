setwd("~/KoborLab/kobor_space/zdong/Population_dataset/MHB_New/Analysis/CMR_analysis")
library(gplots) 
test<-read.table(gzfile("../Old_CMR_analysis/data_CMR.txt.gz"))
test<-test[,1:62]

# http://www.rapidtables.com/web/color/RGB_Color.htm
condition_colors <- unlist(lapply(colnames(test),function(x){
  if(grepl('AFR',x)) '#E69F00' #pink
  else if(grepl('EUR',x)) '#009E73' #grey
  #else if(grepl('ceu',x)) '#008000'
  #else if(grepl('H.s',x)) '#0000FF'
  # else '#FF0000'
}))
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
ggsave(plot = g, width = 2.7, height = 2.7, dpi = 300, filename = "pca_aml.pdf")