require(ggpubr)
require(stringr)
dfm<-read.csv("input",header=T,sep='\t')
dfm$group <- factor(dfm$group, levels = c(1,0))
# dfm$cyl <- factor(dfm$cyl, levels = c('MCC vs SCC','MCC vs NCC'))
p<-ggdotchart(dfm, x = "name", y = "mpg",
              color = "cyl",shape ="group",
              dot.size = 4,
              #palette = "Set1",
              rotate = TRUE,ylab="Fold enrichment",
              sorting = "descending",
              ggtheme = theme_bw()
)+
  # geom_point(aes(shape=factor(group)))+
  scale_shape_manual(values=c(19,1) , guide = "none")+
  theme(legend.position = "top",axis.title.x = element_text(size=12),
        axis.title.y=element_blank(),axis.text = element_text(size=12),
        legend.title =element_blank(),legend.text = element_text(size=12)
  )+
  #scale_y_continuous(breaks = sort(c(seq(0, 20,length.out=5), 1)))+
  scale_x_discrete(labels=function(dfm) str_wrap(dfm,width = 39))+
  scale_y_continuous(limits = c(0, 120))+
  scale_colour_manual(values=c("#e41a1c"))
ggsave("Enrichmen_CpG.pdf",p, width=3.8, height=4)