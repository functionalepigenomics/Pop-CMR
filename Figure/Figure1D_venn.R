

# # Load library
# library(VennDiagram)
# 
# # Generate 3 sets of 200 words
# a<-read.table("../cmr.bed.gz")
# set1 <- a$V6
# b <- read.table("a.txt")
# set2<-unique(b$V6)
# # set3 <- read.table("eur.eas")
# 
# # Prepare a palette of 3 colors with R colorbrewer:
# library(RColorBrewer)
# myCol <-c('#e41a1c','#3182bd')
# 
# # Chart
# venn.diagram(
#   x = list(set1, set2),
#   category.names = c("WGBS","Array"),
#   filename = 'Venn_cor.tif',
#   output=TRUE,
#   
#   # Output features
#   imagetype="tiff" ,
#   height = 800 , 
#   width = 800 , 
#   resolution = 300,
#   compression = "lzw",
#   
#   # Circles
#   lwd = 2,
#   lty = 'blank',
#   fill = myCol,
#   
#   # Numbers
#   cex = .6,
#   fontface = "bold",
#   fontfamily = "sans",
#   
#   # Set names
#   cat.cex = 0.6,
#   cat.fontface = "bold",
#   cat.default.pos = "text",
#   # cat.pos = c(-27, 27, 27),
#   # cat.dist = c(0.055, 0.055, 0.045),
#   cat.fontfamily = "sans",
#   rotation = 1
# )
# 
# draw.pairwise.venn(area1 = length(set1),area2=length(set2))

venn.plot <- draw.pairwise.venn(36657, 416, 416, c("WGBS","Array"),
                                # Circles
                                lwd = 2,
                                lty = 'blank',
                                fill = myCol,
                                # Numbers
                                cex = .6,
                                fontface = "bold",
                                fontfamily = "sans"
                                
);

# Writing to file
tiff(
  filename = tempfile(
    pattern = 'Pairwise_Venn_diagram',
    fileext = '.tiff'
  ),
  compression = "lzw");

grid.draw(venn.plot);
dev.off();