
library(MatrixEQTL)
setwd("~/airr")
args <- commandArgs(trailingOnly=T)
print(args)
## Settings
useModel = modelLINEAR 

SNP_file_name = args[1] 
expression_file_name = args[2] 
if (args[3] == "0")
  covariates_file_name = character() else
  covariates_file_name = args[3] 

output_prefix = args[4] #"hla.9cov.normalize"
prefix = "results/matrixeqtl/"
dir.create(prefix)
output_file_name = paste(prefix, "matrixeqtl.", output_prefix, ".tsv", sep="")
rds_file_name = paste(prefix, "matrixeqtl.", output_prefix, ".rds", sep="") 

print(paste("input file:", SNP_file_name, expression_file_name, covariates_file_name))
print(paste("output file:", rds_file_name, output_file_name))

# Only associations significant at this level will be saved
pvOutputThreshold = as.numeric(args[5])

# Error covariance matrix
errorCovariance = numeric();

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 5000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 5000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Run the analysis
eqtl = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = T);

## Results:
saveRDS(eqtl, file=rds_file_name);
cat('Analysis done in: ', eqtl$time.in.sec, ' seconds', '\n')


