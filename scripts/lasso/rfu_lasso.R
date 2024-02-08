args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
###############################################
### Directories & Variables
# For TEDDY
en.dir <- "results/lasso/"
# For UKB
# en.dir <- "results/lasso_ukbf/"
dir.create(en.dir)
snpset <- args[2] # snpset = "trbsig.ukbflip"
geneset = "tcrnum"

k <- 10 ### k-fold CV

# 1: lasso, 0.5: Elastic Net
alpha <- as.numeric(args[1])

if (length(args) == 2){
  variants.prefix = "data/gwas/sel_"
} else {
  variants.prefix = args[3]
}

library(glmnet)
library(tidyverse)
library(doMC)
registerDoMC(cores = 28)
################################################

whole.expdata.fn = "data/gwas/tcrnum_residuals.tsv"
whole.expdata = t(read.table(whole.expdata.fn, sep="\t", header=T, row.names = 1))
expdata = whole.expdata
snpsets = str_split(snpset, "_")[[1]]
datalist = vector("list", length = length(snpsets))
for (i in 1:length(snpsets)){
  variants.file = paste(variants.prefix, snpsets[i], "_matrixeqtl.traw", sep="")
  datalist[[i]] = as.data.frame(t(read.table(variants.file, sep="\t", header=T, row.names = 1)))
}
X = data.matrix(bind_cols(datalist))
cat("dim of X, ", dim(X), '\n')
print("X")
print(X[1:5,1:5])
print("expdata")
print(expdata[1:5, 1:5])

resultsarray <- array(0,c(dim(expdata)[2],8))
dimnames(resultsarray)[[1]] <- colnames(expdata)
resultscol <- c("gene","alpha","cvm","lambda.iteration","lambda.min","n.snps","R2","pval")
dimnames(resultsarray)[[2]] <- resultscol
workingbest <- en.dir %&% "lasso_alpha" %&% alpha %&% "_" %&% 
  snpset %&% "_" %&% geneset %&% ".txt"
print(workingbest)
if (file.exists(workingbest)){
  print("File exists. Exiting.")
  return()
}
write(resultscol,file=workingbest,ncolumns=8,sep="\t")

weightcol = c("gene","SNP","beta")
workingweight <- en.dir %&% "weights_alpha" %&% alpha %&% "_" %&% snpset %&% 
  "_" %&% geneset %&% ".txt"
print(workingweight)
if (file.exists(workingweight)){
  print("File exists. Exiting.")
  return()
}
write(weightcol,file=workingweight,ncol=3,sep="\t")

set.seed(1001)

for(gene in colnames(expdata)){
  cat(gene %&% '\n')
  exppheno <- expdata[,gene] ### pull expression data for gene

  ##run Cross-Validation over alphalist
  fit <- cv.glmnet(X,exppheno,nfolds=k,alpha=alpha, keep = T, parallel=T) # keep: For R2 calculation
  
  fit.df <- data.frame(fit$cvm,fit$lambda,1:length(fit$cvm)) 
  best.lam <- fit.df[which.min(fit.df[,1]),] 
  cvm.best = best.lam[,1]
  lambda.best = best.lam[,2]
  nrow.best = best.lam[,3] 
  
  ret <- as.data.frame(fit$glmnet.fit$beta[,nrow.best])
  ret[ret == 0.0] <- NA
  bestbetas = as.vector(ret[which(!is.na(ret)),])
  names(bestbetas) = rownames(ret)[which(!is.na(ret))]

  pred.mat <- fit$fit.preval[,nrow.best] 


  if(length(bestbetas) > 0){
    res <- summary(lm(exppheno~pred.mat))
    genename <- as.character(gene)
    rsq <- res$r.squared
    pval <- res$coef[2,4]

    resultsarray[gene,] <- c(genename, alpha, cvm.best, nrow.best, lambda.best, length(bestbetas), rsq, pval)

    bestbetalist <- names(bestbetas)
    betatable<-as.matrix(cbind(bestbetalist,bestbetas))
    betafile<-cbind(genename,betatable[,1],betatable[,2]) ##output "gene","SNP","beta"
    write(t(betafile),file=workingweight,ncolumns=3,append=T,sep="\t") # t() necessary for correct output from write() function

  }else{
    genename <- as.character(gene)
    resultsarray[gene,1] <- genename
    resultsarray[gene,2:8] <- c(NA,NA,NA,NA,0,NA,NA)
  }
  write(resultsarray[gene,],file=workingbest,ncolumns=8,append=T,sep="\t")
}

