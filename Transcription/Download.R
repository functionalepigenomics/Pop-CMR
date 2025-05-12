setwd("~/KoborLab/kobor_space/zdong/Population_dataset/MHB_May3/WGBS_transcript")
# workspace setup
require(oligo)
require(hgu133plus2.db)

### ==== GSE49628
gseID <- 'GSE49628'

# read in and prepare data
# raw CEL files
message('--loading CEL files for ', gseID)

raw <- read.celfiles(
  filenames = list.files('GSE49628_RAW/', pattern = '*cel', full.names = TRUE),
  sampleNames = gsub('\\.cel$', '', list.files('GSE49628_RAW/', pattern = '*cel')))

# RMA
message('--RMA normalising...')
gset <- rma(raw)
message('--Done.')

gset <- exprs(gset)

## 
transcript <- data.frame(gset)
row_NA <- as.numeric(rowSums(is.na(transcript)))
col_NA <- as.numeric(colSums(is.na(transcript)))
transcript$row_NA_ratio <- row_NA/(nrow(transcript)) *100

## Set the cut-off value for each methylation site (now NA portion fewer than 5% are selected) ##56943:10%;56721:5%
transcript_2<-transcript[transcript$row_NA_ratio < 5,]
transcript_2<-transcript_2[1:length(colnames(transcript_2))-1]

## Set the cut-off value for each sample (now NA portion fewer than 1% are selected)
transcript_3 <-rbind(transcript_2,col_NA/length(transcript_2[,1]) *100)
transcript_4<-transcript_3[,transcript_3[length(rownames(transcript_3)),1:(ncol(transcript_3))] < 1]
#transcript_4<-data.frame(TargetID,transcript_4)
transcript_5<-transcript_4[1:nrow(transcript_4)-1,]

gset<-transcript_5
probes <- rownames(gset)
samIDs <- colnames(gset)
# annotate
annotLookup <- select(hgu133plus2.db, keys = probes,
                      columns = c('PROBEID', 'ENSEMBL', 'SYMBOL'))

# remove probes with any NA mapping
annotLookup <- annotLookup[!is.na(annotLookup$ENSEMBL) & !is.na(annotLookup$SYMBOL),]
annotLookup <- annotLookup[!duplicated(annotLookup$PROBEID),]

# look up the ensembl ID and gene symbol
probes <- probes[which(probes %in% annotLookup$PROBEID)]
gset <- gset[probes,]
all(rownames(gset)==probes)
all(probes == annotLookup[match(probes, annotLookup$PROBEID),'PROBEID'])
geneid <- annotLookup[match(probes, annotLookup$PROBEID),'SYMBOL']
ens <- annotLookup[match(probes, annotLookup$PROBEID),'ENSEMBL']

# finalise the dataset
final <- data.frame(ens = ens, symbol = geneid, gset)
head(final)

### ==== GSE11582
library(GEOquery)
getGEOSuppFiles("GSE11582")
setwd("GSE11582/")
untar("GSE11582_RAW.tar")
setwd("~/KoborLab/kobor_space/zdong/Population_dataset/MHB_May3/WGBS_transcript")

gseID_1<- 'GSE11582'

# read in and prepare data
# raw CEL files
message('--loading CEL files for ', gseID_1)

raw_1 <- read.celfiles(
  filenames = list.files('GSE11582/', pattern = '*CEL', full.names = TRUE),
  sampleNames = gsub('\\.CEL$', '', list.files('GSE11582/', pattern = '*CEL'))
  )

# RMA
message('--RMA normalising...')
gset_1 <- rma(raw_1)
message('--Done.')

gset_1 <- exprs(gset_1)

transcript <- data.frame(gset_1)
row_NA <- as.numeric(rowSums(is.na(transcript)))
col_NA <- as.numeric(colSums(is.na(transcript)))
transcript$row_NA_ratio <- row_NA/(nrow(transcript)) *100

## Set the cut-off value for each methylation site (now NA portion fewer than 5% are selected) ##56943:10%;56721:5%
transcript_2<-transcript[transcript$row_NA_ratio < 5,]
transcript_2<-transcript_2[1:length(colnames(transcript_2))-1]

## Set the cut-off value for each sample (now NA portion fewer than 1% are selected)
transcript_3 <-rbind(transcript_2,col_NA/length(transcript_2[,1]) *100)
transcript_4<-transcript_3[,transcript_3[length(rownames(transcript_3)),1:(ncol(transcript_3))] < 1]
#transcript_4<-data.frame(TargetID,transcript_4)
transcript_5<-transcript_4[1:nrow(transcript_4)-1,]

gset_1<-transcript_5

probes_1 <- rownames(gset_1)
samIDs_1 <- colnames(gset_1)

# annotate
annotLookup <- select(hgu133plus2.db, keys = probes_1,
                      columns = c('PROBEID', 'ENSEMBL', 'SYMBOL'))

# remove probes with any NA mapping
annotLookup <- annotLookup[!is.na(annotLookup$ENSEMBL) & !is.na(annotLookup$SYMBOL),]
annotLookup <- annotLookup[!duplicated(annotLookup$PROBEID),]

# look up the ensembl ID and gene symbol
probes_1 <- probes_1[which(probes_1 %in% annotLookup$PROBEID)]
gset_1 <- gset_1[probes_1,]
all(rownames(gset_1)==probes_1)
all(probes_1 == annotLookup[match(probes_1, annotLookup$PROBEID),'PROBEID'])
geneid_1 <- annotLookup[match(probes_1, annotLookup$PROBEID),'SYMBOL']
ens_1 <- annotLookup[match(probes_1, annotLookup$PROBEID),'ENSEMBL']

# finalise the dataset
final_1 <- data.frame(ens = ens_1, symbol = geneid_1, gset_1)
head(final_1)

## Combine
affyA<-gset
affyB<-gset_1
pNListA = rownames(affyA)
pNListB = rownames(affyB)
subsetList = pNListA[pNListA %in% pNListB]

subsetPmA = affyA[unique(subsetList),]
subsetPmB = affyB[unique(subsetList),]
allPm = cbind(subsetPmA, subsetPmB)

# verbose = TRUE
# background = TRUE
# normalize = TRUE
# bgversion = 2
# ngenes = length(unique(subsetList))
# pNList = split(0:(length(subsetList) - 1), subsetList)
# exprs <- .Call("rma_c_complete", allPm, 
#                pNList, ngenes, normalize, background, bgversion, 
#                verbose, PACKAGE = "affy")

gset_2 <- allPm
probes_2 <- rownames(allPm)
samIDs_2 <- colnames(allPm)

# annotate
annotLookup <- select(hgu133plus2.db, keys = probes_2,
                      columns = c('PROBEID', 'ENSEMBL', 'SYMBOL'))

# remove probes with any NA mapping
annotLookup <- annotLookup[!is.na(annotLookup$ENSEMBL) & !is.na(annotLookup$SYMBOL),]
annotLookup <- annotLookup[!duplicated(annotLookup$PROBEID),]

# look up the ensembl ID and gene symbol
probes_2 <- probes_2[which(probes_2 %in% annotLookup$PROBEID)]
gset_2 <- gset_2[probes_2,]
all(rownames(gset_2)==probes_2)
all(probes_2 == annotLookup[match(probes_2, annotLookup$PROBEID),'PROBEID'])
geneid_2 <- annotLookup[match(probes_2, annotLookup$PROBEID),'SYMBOL']
ens_2 <- annotLookup[match(probes_2, annotLookup$PROBEID),'ENSEMBL']

# finalise the dataset
colnames(gset_2)<-samIDs_2
final_2 <- data.frame(ens = ens_2, symbol = geneid_2, gset_2)
head(final_2)
write.table(final_2,file="Combined_needadjust_batch.txt",quote=F,sep='\t',row.names = T,col.names = T)

