setwd("~/KoborLab/kobor_space/zdong/Population_dataset/MHB_New/Analysis/twoloci_cor/Genomicfeatures/Density/")
require(ggpubr)
dfm<-read.csv("input",header=T,sep="\t")
# p<-ggdotchart(dfm, x = "name", y = "mpg",
#               #group = "cyl", 
#               color = "cyl",
#               dot.size = 4,
#               palette = "Set1",
#               rotate = TRUE,ylab="Ratio of Genetic-associated CpGs across diverse CpG types",
#               sorting = "descending",
#               ggtheme = theme_bw()
# )+
#   theme(legend.position = c(0.88,0.15),
#         axis.title.y=element_blank(),
#         legend.text = element_text(size=12),
#         axis.title.x = element_text(size=12),
#         axis.text = element_text(size=11),
#         legend.title =element_blank()
#   )+#geom_hline(yintercept=0, linetype="dashed", color = "black",show.legend=F)+
#   scale_y_continuous(breaks = sort(c(seq(0, 0.35,length.out=8), 1)))
# ggsave("lollipop.chart.pdf",p, width=2, height=1.3, units="in", scale=3)

dfm$cyl <- factor(dfm$cyl, levels = dfm$cyl)
p<-ggplot(data=dfm, aes(x=cyl,y=mpg,fill=name))  +theme_classic2()+#scale_y_continuous(limits = c(0,0.6),breaks = c(0.1,0.2,0.3,0.4,0.5,0.6))+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12), 
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=12),
        axis.text = element_text(size=12),
        legend.text=element_blank(),
        legend.title =element_blank()
  )+ylab("Fold enrichment")+
  #geom_rect(aes(),xmin =-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,alpha = 0.1) +
  geom_bar(aes(fill = name),stat = "identity",position = "dodge", width=0.8) +
  scale_fill_manual(values=c("#e41a1c","#43a2ca","grey"))+ 
  coord_flip()+scale_y_continuous(limits = c(0, 8))+
  geom_hline(yintercept = 1,linetype = "dashed")
ggsave("CMR_genomicfeatures.pdf",p, width=5, height=5.4)#, units="in", scale=3)

