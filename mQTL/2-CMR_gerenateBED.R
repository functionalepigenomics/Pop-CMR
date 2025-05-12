eur <-read.table('../../Array_methylation/MHB_recall/AML/Genomericfeatures/pop-MHB1.sorted.txt')
eur$geneid<-paste(eur$V1,eur$V2,eur$V3,sep = ':')
colnames(eur)[1:3]<-c('chr','s1','s2')
eur<-eur[,c('geneid','chr','s1','s2')]
eur$chr<-paste('chr',eur$chr,sep='')
write.table(eur,file='MHB_coordinate1.bed',quote = F,sep='\t',row.names = F)
