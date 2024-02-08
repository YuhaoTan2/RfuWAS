setwd("~/airr")
"%&%" = function(a,b) paste0(a,b)
library(tidyverse)
library(glmnet)
library(vroom)
library(doParallel)
library(foreach)
args <- commandArgs(trailingOnly=T)

# Set the number of cores to be used
num_cores <- 2

# Set up parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# For TEDDY
# en.dir <- "results/lasso/"
# For UKB
en.dir <- "results/lasso_ukbf/"
snplabel = ".ukbflip" # ".ukb", ".tesla", ".intersect"
model.dir = "models/lasso"%&% snplabel %&% "/"
test.set = "ukb_hg19"

all.model.best = read.table(en.dir %&% "bestmodel"%&% snplabel %&%".txt",header=T,sep="\t")
alpha <- 1
associated.genes = all.model.best$gene
cat("Predicting", length(associated.genes), "genes\n")

variants.file = "data/ukb/lasso_variant_sets/ukb_sampleqc_variantqc_matrixeqtl.traw"

# Split variants.file to parts
num_parts <- 12 * num_cores
num_cols <- ncol(vroom(variants.file, n_max = 1))
cat("Total sample number", num_cols, '\n')
cols_per_part <- ceiling(num_cols / num_parts)

predict_part <- function(i) {
  log_fn = paste0("log/predict_ukb_",i,".out")
  cat <- function(...) base::cat(..., file = log_fn, append = TRUE)
  sink(log_fn, append = TRUE)
  start_col <- max(2, (i - 1) * cols_per_part + 1)
  end_col <- min(i * cols_per_part, num_cols)
  cat("Calculating part", i, start_col, end_col, '\n')

  datalist = vroom(variants.file, col_select = start_col:end_col)
  X.ref = data.matrix(t(datalist))
  cat(dim(X.ref), '\n')

  # write sample name in predict file
  workingpred = en.dir %&% "bestmodel"%&% snplabel %&%"_"%&% test.set %&%"_"%&% i %&%".txt"
  print(workingpred)
  if (file.exists(workingpred)) {
  message("The file exists.")
  return()
  }
  pred.cols = nrow(X.ref) + 1
  X.ref.IID = gsub(".*_","",rownames(X.ref))
  write(c("gene", X.ref.IID),file=workingpred,ncolumns=pred.cols,sep="\t")

  for(gene in associated.genes){
    cat(gene %&% '\n')
    fit = readRDS(model.dir %&% "bestmodel_" %&% gene %&%".rds")
    pred.mat= predict(fit, newx = X.ref)
    write(c(as.character(gene), t(pred.mat)),file=workingpred,ncolumns=pred.cols,append=T,sep="\t")
  }
}

# Get the list of all objects in the global environment
all_objects <- ls(globalenv())

# Export all objects to the cluster
clusterExport(cl, all_objects)

clusterEvalQ(cl, lapply(c("tidyverse", "glmnet",
                "vroom"), library, character.only = TRUE))

# Use foreach to execute the loop in parallel
foreach(i = (num_cores * (as.numeric(args[1]) - 1) + 1) : (num_cores * as.numeric(args[1]))) %dopar% {
  predict_part(i)
}

# Stop the parallel backend
stopCluster(cl)
registerDoSEQ()

