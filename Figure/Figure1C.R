library(corrplot)
library(RColorBrewer)
test<-read.table("../MHB.CpG1.txt")
mtcars<-test[which(test$V1 %in% c('10:16479126:16479166')),c(2:ncol(test))]
# mtcars<-mtcars[1:14,]
rownames(mtcars)<-mtcars$V2
mtcars<-mtcars[,2:ncol(mtcars)]
# mtcars<-mtcars[,colSums(is.na(mtcars)) != nrow(mtcars),]
mtcars<-mtcars[ , colSums(is.na(mtcars))==0]
M <-cor(t(mtcars))


pdf(file = "cc.pdf",width=6, height=6)
cex.before <- par("cex")
par(cex = 0.7)
cl<- c("#7C8489","#4fB3A4","#F5B977") # "#56A36C","#5E8579","#2E68AA",
col2 <- colorRampPalette(cl)
#col2<- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M, outline = T, addgrid.col = "white",mar=c(2,2,2,2),col = col2(10),
         cl.pos = "r", tl.col = "indianred4", tl.srt=45,#method ="ellipse",
         tl.cex = 1/par("cex"),addCoef.col = "black",col.lim = c(0,1),is.corr = F,
         cl.cex = 1/par("cex"), #addCoefasPercent = TRUE,
         type = "upper", tl.pos = "tl", bg="azure2") #azure2
dev.off()