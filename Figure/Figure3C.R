meth<-read.table("../../DM-CMR.txt.gz")
target<-read.table("../mQTL_output_lead.txt.gz",header = T)
meth1<-meth[target$gene,]
meth2<-meth[-which(rownames(meth) %in% target$gene),]
summary(abs(meth1$dif))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1648  0.3210  0.3666  0.3886  0.4459  0.7821
summary(abs(meth2$dif))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1377  0.2899  0.3268  0.3231  0.3622  0.4925 
wilcox.test(abs(meth1$dif),abs(meth2$dif))$p.value # 0.002089807

## ====== boxplot
t1<-abs(meth1$dif)
t1<-data.frame(t1,"SNP-associated Pop-CMRs")#"Population-specific effect and allelic frequency")
colnames(t1)<-c("AML changes",'Type')
t2<-abs(meth2$dif)
t2<-data.frame(t2,'Other Pop-CMRs')#"Population-specific effect not allelic frequency")
colnames(t2)<-c("AML changes",'Type')
data<-rbind(t2,t1)
colnames(data)<-c("AML_changes",'Type')
data$Type <- factor(data$Type, levels = c("SNP-associated Pop-CMRs", "Other Pop-CMRs"))

library(readr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(reshape2)
raincloud_theme = theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  #axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.title=element_text(size=16),
  legend.text=element_text(size=16),
  legend.position = "right",
  plot.title = element_text(lineheight=.8, face="bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
g <- ggplot(data, aes(y = AML_changes, x = Type, fill =Type)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = AML_changes, color = Type), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 3) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  #scale_color_brewer(palette = "Set1") +
  #scale_fill_brewer(palette = "Set1") +
  scale_fill_manual(values = c("#e41a1c","#fc9272"))+
  scale_color_manual(values = c("#e41a1c","#fc9272"))+
  # coord_flip() +
  theme_bw() +
  raincloud_theme+theme(panel.grid = element_blank(),legend.position = "none", 
                        panel.background = element_rect(fill = 'transparent', color = 'black'), 
                        legend.title = element_blank(), legend.key = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        strip.text = element_text(size=14), 
                        axis.title.x=element_blank(),
                        axis.title.y=element_text(size=14),
                        axis.text = element_text(size=14))+ylim(0, 1)
ggsave("AML_comparsion.pdf",g, width=4, height=2.6)#, units="in", scale=3)
