num<-read.table(gzfile("../../twoloci_cor/Permutation/Permutation/pair_merged.bed.gz"))
num$V4=num$V4+1
num<-num[num$V4>2,]
rownames(num)<-paste(num$V1,num$V2,num$V3,sep=":")
pop<-read.table("../DM-CMR.txt")
summary(num$V4)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.000   3.000   3.000   3.667   4.000 146.000 
a<-num[rownames(pop),]
summary(a$V4)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.000   3.000   3.000   3.891   4.000  16.000


n=0
for (i in 1:1000){
  b<-num$V4[sample(length(num$V4),length(a$V4))]
  if(mean(b)>=mean(a$V4)){
    n=n+1
  }
}

n/1000 # 0.095

num$name<-"CMRs"
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
g <- ggplot(num, aes(x=name,y = V4,fill=name)) +
  # geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  # geom_point(aes(y = AML_changes, color = Type), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = .6, guides = FALSE, outlier.shape = NA, alpha = 1) +
  stat_summary(fun.y = mean, geom = "point", shape=20, size=3, color="black", fill="black",
               position = position_dodge(1.1))+
  # expand_limits(x = 3) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  #scale_color_brewer(palette = "Set1") +
  #scale_fill_brewer(palette = "Set1") +
  scale_fill_manual(values = c("#e41a1c"))+
  scale_color_manual(values = c("#e41a1c"))+
  coord_flip() +
  theme_bw() +ylab("CpG count")+xlab("")+
  raincloud_theme+theme(panel.grid = element_blank(),legend.position = "none", 
                        panel.background = element_rect(fill = 'transparent', color = 'black'), 
                        legend.title = element_blank(), legend.key = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        strip.text = element_text(size=14), 
                        axis.title.x=element_text(size=14),
                        axis.title.y=element_blank(),
                        axis.text = element_text(size=14))+coord_flip(ylim = c(0, 10))
ggsave("count.pdf",g, width=4, height=1)#, units="in", scale=3)

##### ====================== LENGTH ==========================================
### length
num$length<-num$V3-num$V2+1
a<-num[rownames(pop),]
summary(num$length)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 5.00   16.00   28.00   50.81   55.00 3934.00
summary(a$length)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.0    25.0    58.0   109.5   139.0   852.0 

n=0
for (i in 1:1000){
  b<-num$length[sample(length(num$length),length(a$length))]
  if(mean(b)>=mean(a$length)){
    n=n+1
  }
}

n/1000 # 0

g <- ggplot(num, aes(x=name,y = length,fill=name)) +
  # geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  # geom_point(aes(y = AML_changes, color = Type), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = .6, guides = FALSE, outlier.shape = NA, alpha = 1) +
  stat_summary(fun.y = mean, geom = "point", shape=20, size=3, color="black", fill="black",
               position = position_dodge(1.1))+
  # expand_limits(x = 3) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  #scale_color_brewer(palette = "Set1") +
  #scale_fill_brewer(palette = "Set1") +
  scale_fill_manual(values = c("#e41a1c"))+
  scale_color_manual(values = c("#e41a1c"))+
  theme_bw() +ylab("Length")+xlab("")+
  raincloud_theme+theme(panel.grid = element_blank(),legend.position = "none", 
                        panel.background = element_rect(fill = 'transparent', color = 'black'), 
                        legend.title = element_blank(), legend.key = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        strip.text = element_text(size=14), 
                        axis.title.x=element_text(size=14),
                        axis.title.y=element_blank(),
                        axis.text = element_text(size=14))+coord_flip(ylim = c(0, 100))
ggsave("length.pdf",g, width=4, height=1)#, units="in", scale=3)

#================================== bar chart ============================================
library(ggplot2)
dfm<-data.frame(name=c("CMRs"),mpg=c(36657),cyl=c("CMRs"))
# dfm$name <- factor(dfm$name, levels = c("SNPs for Pop-CMRs","Random SNPs"))
p<-ggplot(data=dfm, aes(x=name,y=mpg))  +theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=14),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text = element_text(size=14),
        legend.title =element_blank()
  )+ylab("Number")+xlab("")+
  #geom_rect(aes(),xmin =-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,alpha = 0.1) +
  geom_bar(aes(fill = name),stat = "identity",position = "dodge") + ylim(0,40000)+
  scale_fill_manual(values=c("#e41a1c"))+raincloud_theme+theme(panel.grid = element_blank(),legend.position = "none", 
                                                               panel.background = element_rect(fill = 'transparent', color = 'black'), 
                                                               legend.title = element_blank(), legend.key = element_blank(),
                                                               panel.grid.major = element_blank(),
                                                               panel.grid.minor = element_blank(),
                                                               strip.text = element_text(size=14), 
                                                               axis.title.x=element_text(size=14),
                                                               axis.title.y=element_blank(),
                                                               axis.text = element_text(size=14))+
  coord_flip()
#facet_wrap(~ tissue, nrow=2)
ggsave("Number.pdf",p, width=4, height=1, dpi = 300)
