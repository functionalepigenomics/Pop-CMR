a<-read.table("pheno.txt",row.names = 1,header=T)
rownames(a)<-sub("[_].*", "", rownames(a))
rownames(a)<-sub("[.].*", "", rownames(a))
b<-read.table("README.txt",row.names = 1,header=T)
a$Sample<-NA
sum(is.na(a$Sample)) # 305
for (i in 1:nrow(a)){
  for (j in 1: nrow(b)){
    if (rownames(a[i,]) == as.character(b$Tans_ID[j])){a$Sample[i]<-rownames(b[j,])}
  }
}
sum(is.na(a$Sample)) # 196

b<-read.table("../Array_transcript/a.log",row.names = 1,header=T,sep='\t')
b<-data.frame(t(b))
b$Sample_id<-sub('GM', 'NA', b$Sample_id)
sum(is.na(a$Sample)) # 196
for (i in 1:nrow(a)){
  for (j in 1: nrow(b)){
    if (rownames(a[i,]) == rownames(b[j,])){a$Sample[i]<-as.character(b$Sample_id[j])}
  }
}
sum(is.na(a$Sample)) # 5
a$Sample[is.na(a$Sample)]<-c("B1","B2","B3","NA12878","NA10860")
sum(is.na(a$Sample)) # 0
length(unique(a$Sample))==nrow(a)
# [1] TRUE

sex<-read.table("../Array_methylation/Meth_array_sampleinfo.txt",header=T,row.names = 1)
a$sex<-NA
sum(is.na(a$sex)) # 305
for (i in 1:nrow(a)){
  for (j in 1:nrow(sex)){
    if (as.character(a$Sample[i]) == as.character(sex$Sample[j])){a$sex[i]<-as.numeric(sex$sex[j])}
  }
}
sum(is.na(a$sex)) # 10
a$sex[is.na(a$sex)][4:10]<-c(1,2,2,2,2,2,1)
a$Sample[is.na(a$sex)]
# [1] "B1" "B2" "B3"
write.table(a,file="Transcript_sampleinfo.txt",quote = T,sep='\t')

