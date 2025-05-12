# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# 
# Be sure to use an up to date version of R and Matrix eQTL.

#setwd("~/KoborLab/kobor_space/kandy/home/zdong/Population_dataset/MHBmap/SNP_new_forimputation/Correation/AFR")
library(MatrixEQTL)

## Location of the package with the data files.
#base.dir = find.package('MatrixEQTL');
base.dir = '~/KoborLab/kobor_space/zdong/Population_dataset/MHB_May3/Genotype/BackWGBS_AFR';

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste(base.dir, "/Input_genotype_allchr.txt", sep="");
snps_location_file_name = paste(base.dir, "/SNP_coordinate1.bed", sep="");

# Gene expression file name
expression_file_name = paste(base.dir, "/Input_methylation.txt", sep="");
gene_location_file_name = paste(base.dir, "/MHB_coordinate1.bed", sep="");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = character() # paste(base.dir, "/Input_covariate.txt", sep="");

# Output file name
output_file_name_cis = "result_cis.txt";
output_file_name_tra = "result_tra.txt";

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 0.05;
pvOutputThreshold_tra = 0.05;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e5;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_cis,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

# unlink(output_file_name_tra);
# unlink(output_file_name_cis);

## Results:

# cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
# cat('Detected local eQTLs:', '\n');
# show(me$cis$eqtls)
# cat('Detected distant eQTLs:', '\n');
# show(me$trans$eqtls)

## Plot the Q-Q plot of local and distant p-values

plot(me)










