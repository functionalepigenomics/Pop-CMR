#### 2. Run MatrixEQTL Analysis ####

# Location of the data files.
base.dir = "~/KoborLab/kobor_space/zdong/Population_dataset/MHB_New/Analysis/CMR_analysis/Genotype"; # Replace with your working directory

## Settings

# Genotype and variant information file names
SNP_file_name = paste(base.dir,"/BSGS_Genotype.txt",sep="");
snps_location_file_name = paste(base.dir,"/SNP_Pos.txt",sep="");

snps <- read.table(SNP_file_name, sep="\t", header=T, row.names = 1, check.names = F)
snps[1,]

# Methylation and probe information file names
expression_file_name = paste(base.dir,"/BSGS_Methylation.txt",sep="");
gene_location_file_name = paste(base.dir,"/Probe_Pos.txt",sep="");

expr = read.table(expression_file_name, sep="\t", header=T, row.names = 1, check.names = F)
expr[1,]

# Covariate file name
covariates_file_name = paste(base.dir,"/BSGS_covariates.txt",sep="");
cvrt = read.table(covariates_file_name, sep = "\t", header=T, row.names = 1, check.names = F)
cvrt[1,]

#### 3. Plotting mQTLs ####

probe<-expr["6:159703533:159703779",]
genotypes<-snps["6:159281319:A:G",]
cpg<-rownames(probe)
variant<-"rs543547"
# Merge genotype and methylation data
data <- data.frame(t(genotypes),t(probe))
names(data)<-c("genotype","probe")
head(data)
data$genotype[data$genotype==0]<-"AA"
data$genotype[data$genotype==1]<-"AG"
data$genotype[data$genotype==2]<-"GG"

# Plot
library(ggplot2)
plot <- ggplot(data, aes(x = genotype, y = probe)) + geom_boxplot(fill = "#E41A1C") + ylab("AML") + xlab("Genotype") + theme(text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) + geom_point()
plot
pdf(paste(cpg,"_",variant,"_plot.pdf",sep=""),width = 3,height = 2.4)
plot
dev.off()