test<-read.table("Permutation/MHB.CpG1.txt.gz")
# mtcars<-test[,c(2:ncol(test))]
# 
# rownames(mtcars)<-mtcars$V2
# mtcars<-mtcars[,2:ncol(mtcars)]
# # mtcars<-mtcars[,colSums(is.na(mtcars)) != nrow(mtcars),]
# mtcars<-mtcars[ , colSums(is.na(mtcars))==0]

p<-c(0,0);cor<-c(0,0);n=0;block<-c('1','1')
for (i in 1:(nrow(test)-1)) {
  if (as.character(test$V1[i]) == as.character(test$V1[i+1])){n=n+1;
  cal<-cor.test(as.numeric(test[i,c(3:64)]),as.numeric(test[(i+1),c(3:64)]),na.action = "na.exclude", alternative = c("two.sided"),
                method = c("pearson"))
  block[n]<-as.character(test$V1[i])
  cor[n]<-cal$estimate
  p[n]<-cal$p.value
  }
}
afr<-data.frame(block,cor,p)
write.table(afr,"recalculated_cor.txt",quote = F,sep="\t")



a<-read.table(gzfile("all.txt.gz"),header=T)
# colnames(a)<-c("block","cor","p","fdr")
nrow(a) # 14000010
a[c('chr', 'first','second')] <- str_split_fixed(a$block, ':', 3)
a$length<-as.numeric(a$second)-as.numeric(a$first)+1
a<-a[a$length<=400,] # 12049591
a$fdr<-p.adjust(a$p,method="fdr") 
pair<-a[a$cor>0.5 & a$fdr<0.05,] # 726893
write.table(pair,file="pair_withcor.txt",sep="\t",quote = F)

pair<-read.table("pair_withcor.txt",header = T)
require(ggplot2)
theme_set(theme_classic())
g <- ggplot(pair,aes(cor))+
  geom_density(fill='#e41a1c',color="#e41a1c",alpha = 0.8#,position = "stack"
  )+ labs(x="Correlation coefficient r")+
  # geom_rug(aes(x = cor, y = 0,colour=factor(num)),
  #          position = position_jitter(height = 0))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y  = element_blank(),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size=11),
        legend.position=c(0.12,0.9),legend.title =  element_blank(),
        legend.spacing.x = unit(0.1, 'cm'),
        axis.text.y = element_blank(),
        axis.text.x =element_text(size=12),
        axis.ticks.x = element_line(colour = "black", size = 0.8),
        axis.ticks.y=element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())+scale_fill_manual(values=c('#e41a1c'))+
  geom_vline(aes(xintercept=0.5),color="black", linetype="dashed", size=0.8)
# scale_color_manual(values=c("#E69F00","#009E73"))
ggsave(plot = g, width = 3, height = 1.5, dpi = 300, filename = "density_cmr.pdf")