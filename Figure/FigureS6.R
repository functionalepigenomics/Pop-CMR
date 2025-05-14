library(ggVennDiagram)
library(ggplot2)

MSCC<- read.table("afr.eas")
SCC<-read.table("eur.afr")
NCC<-read.table("eur.eas")
MSCC1<- read.table("afr.sas")
SCC1<-read.table("eur.sas")
NCC1<-read.table("eur.amr")
MSCC2<- read.table("afr.amr")
SCC2<-read.table("amr.sas")
NCC2<-read.table("eas.sas")
MSCC3<-read.table("eas.amr")

overlap<-list(AFRvsEAS=MSCC$V1,EURvsAFR=SCC$V1,
              EURvsEAS=NCC$V1,
              AFRvsSAS=MSCC1$V1,EURvsSAS=SCC1$V1,
              EURvsAMR=NCC1$V1,
              AFRvsAMR=MSCC2$V1,AMRvsSAS=SCC2$V1,
              EASvsSAS=NCC2$V1,EASvsAMR=MSCC3$V1
              )
library("UpSetR")
library("ComplexHeatmap")
m = make_comb_mat(overlap)
pdf(file = "upset_fivepopulations.pdf",width = 9, height = 4)
UpSet(m, top_annotation = upset_top_annotation(m, 
                                               annotation_name_rot = 90,
                                               annotation_name_side = "right",
                                               axis_param = list(side = "right")))
dev.off()