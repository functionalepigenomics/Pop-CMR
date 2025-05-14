library(ggVennDiagram)
library(ggplot2)

MSCC<- read.table("afr.eas")
SCC<-read.table("eur.afr")
NCC<-read.table("eur.eas")


overlap<-list(AFRvsEAS=MSCC$V1,EURvsAFR=SCC$V1,
              EURvsEAS=NCC$V1)
library("UpSetR")
library("ComplexHeatmap")
m = make_comb_mat(overlap)
pdf(file = "upset.pdf",width = 4, height = 2.2)
UpSet(m, top_annotation = upset_top_annotation(m, 
                                               annotation_name_rot = 90,
                                               annotation_name_side = "right",
                                               axis_param = list(side = "right")))
dev.off()