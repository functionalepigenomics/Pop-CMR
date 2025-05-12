#================================== bar chart ============================================
library(ggplot2)
dfm<-data.frame(name=c("Pop-CMRs","Random CMRs"),mpg=c(68.3,57.5),cyl=c("Pop-CMRs","Random CMRs"))
dfm$name <- factor(dfm$name, levels = c("Pop-CMRs","Random CMRs"))
p<-ggplot(data=dfm, aes(x=name,y=mpg))  +theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        axis.text = element_text(size=12),
        legend.title =element_blank()
  )+ylab("Percentage of high\nmethylation in AFR (%)")+xlab("")+
  #geom_rect(aes(),xmin =-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,alpha = 0.1) +
  geom_bar(aes(fill = name),stat = "identity",position = "dodge") + ylim(0,100)+
  scale_fill_manual(values=c("#e41a1c","#636363"))+
  theme(axis.text.x = element_text(angle = 20))
#facet_wrap(~ tissue, nrow=2)
ggsave("percentage of hyper-methyaltion.pdf",p, width=2, height=3, dpi = 300)