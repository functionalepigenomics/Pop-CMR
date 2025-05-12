setwd("~/KoborLab/kobor_space/zdong/Population_dataset/MHB_New/Analysis/CMR_analysis/Genomicfeatures")
a<-read.csv('cmr_pair.count.txt',sep = '\t',header=F)

files <- list.files(path="./", pattern="*.cmr_background.count.txt", full.names=TRUE, recursive=FALSE)
i=0;fold<-c(1,1);p<-c(1,1);
for (x in 1:length(files)) {i=i+1;
b<-read.table(files[x],header=F)
fold[i]<-a$V1[i]/mean(b$V1)
fold #
c<-b$V1[b$V1 >= a$V1[i]]
p[i]<-length(c)/1000000 #
}
fd<-read.table('nohup.out')
total<-cbind(as.character(fd$V1),a$V1,fold,p)
colnames(total)[1:2]<-c('Feature','Num.')
total
write.table(total,file='Pop-MHB_enrichment.txt',quote=F,sep='\t',row.names=F)


#### ==== deleption
files <- list.files(path="./", pattern="*.cmr_background.count.txt", full.names=TRUE, recursive=FALSE)
i=0;fold<-c(1,1);p<-c(1,1);
for (x in 1:length(files)) {i=i+1;
b<-read.table(files[x],header=F)
fold[i]<-a$V1[i]/mean(b$V1)
fold #
c<-b$V1[b$V1 <= a$V1[i]]
p[i]<-length(c)/2000 #
}
fd<-read.table('nohup.out')
total<-cbind(as.character(fd$V1),a$V1,fold,p)
colnames(total)[1:2]<-c('Feature','Num.')
total


# Feature                               Num. fold                p       
# [1,] "10_Txn_Elongation.bed"               "3"  "1.16981867810489"  "0.748" 
# [2,] "1_Active_Promoter.bed"               "4"  "0.446229361892012" "0.048" 
# [3,] "2_Weak_Promoter.bed"                 "4"  "0.59889204970804"  "0.192" 
# [4,] "3_Poised_Promoter.bed"               "0"  "0"                 "0.0955"
# [5,] "8_Insulator.bed"                     "4"  "2.02685583987839"  "0.955" 
# [6,] "9_Txn_Transition.bed"                "2"  "1.03761348897536"  "0.692" 
# [7,] "cpgIslandExt.hg19.bed"               "4"  "0.35885704032656"  "0.011" 
# [8,] "gencode.v19.2wayconspseudos.bed"     "1"  "1.50829562594268"  "0.8565"
# [9,] "gencode.v19.long_noncoding_RNAs.bed" "13" "1.16555341372663"  "0.779" 
# [10,] "human.refgene.3utr.bed"              "4"  "1.38600138600139"  "0.828" 
# [11,] "human.refgene.5utr.bed"              "4"  "1.42500890630566"  "0.8455"
# [12,] "human.refgene.codingexon.bed"        "5"  "0.66894106629206"  "0.225" 
# [13,] "human.refgene.downstream1kb.bed"     "3"  "1.50451354062187"  "0.8525"
# [14,] "human.refgene.intron.bed"            "58" "1.19319467588307"  "0.9705"
# [15,] "Strong_Enhancer.bed"                 "11" "0.641305932079872" "0.064" 
# [16,] "test.real.shelves.bed"               "4"  "0.810619110345526" "0.453" 
# [17,] "test.real.shores.bed"                "10" "0.605895361871005" "0.037" 
# [18,] "Weak_Enhancer.bed"                   "19" "1.65339598833921"  "0.988" 