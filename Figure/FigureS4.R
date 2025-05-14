library(gplots)
library(ggpubr)
probenames<-re_w$re.V1.re.p...0.05.
gene<-c("GAS7","SASH1")
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
