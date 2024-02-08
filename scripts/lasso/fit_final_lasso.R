# Fit final lasso model after cross-validation
# Models saved to "models/lasso.ukb/bestmodel_" %&% gene %&%".rds"
# Weights saved to "results/lasso/weights_bestmodel.ukb.txt"
setwd("~/airr")
"%&%" = function(a,b) paste0(a,b)
library(tidyverse)
library(readr)
library(glmnet)

# For TEDDY
# en.dir <- "results/lasso/"
# For UKB
en.dir <- "results/lasso_ukbf/"
snplabel = ".ukbflip" # ".ukb", ".tesla", ".intersect" ".emerson" 
model.dir = "models/lasso"%&% snplabel %&% "/"
variants.prefix = "data/ukb/lasso_variant_sets/"
if (snplabel %in% c(".ukb", ".ukbflip", ".emerson")){
  ref.snpset.name = paste0("trb.hla", snplabel)
}  else {stop("ref.snpset.name not defined. 
             Please carefully check the parameter, 
             or you may get wrong order of variants!")}
#  ref.snpset.name = paste0("trb", snplabel, "_hla", snplabel)

snplabel.for.variantfile = snplabel
dir.create(model.dir, recursive = TRUE)

all.model.best = read.table(en.dir %&% "bestmodel"%&% snplabel %&%".txt",header=T,sep="\t")
cat("RFUs with lasso model", dim(all.model.best)[1], "\n")
alpha <- 1
associated.genes = all.model.best$gene
# 659 samples * 1205 genes 
whole.expdata.fn = "data/gwas/tcrnum_residuals.tsv"
whole.expdata = t(read.table(whole.expdata.fn, sep="\t", header=T, row.names = 1))
expdata = whole.expdata[, as.character(associated.genes)]
print(dim(expdata))

all.snpsets = unique(all.model.best$modelname)
all.model.best$modelname = factor(all.model.best$modelname, all.snpsets, 
                                 seq_along(all.snpsets))
all.model.best <- data.frame(all.model.best[,-1], row.names = all.model.best[,1])

get.X.ref = function(){
  snpsets = str_split(ref.snpset.name, "_")[[1]]
  datalist = vector("list", length = length(snpsets))
  for (i in 1:length(snpsets)){
    variants.file = paste(variants.prefix, snpsets[i], "_matrixeqtl.traw", sep="")
    datalist[[i]] = as.data.frame(t(read.table(variants.file, sep="\t", header=T, row.names = 1)))
  }
  return(data.matrix(bind_cols(datalist)))
}
X.ref = get.X.ref()
print("X")
print(X.ref[1:5,1:5])
print("expdata")
print(expdata[1:5, 1:5])

all.X.exclude = vector("list", length = length(all.snpsets))
for (j in seq_along((all.snpsets))){
  # 658 samples * 4323 variants
  if (all.snpsets[j] == "trb_hla" & snplabel == ".ukb") {
    snpsets = "trb.hla"
  } else snpsets = str_split(all.snpsets[j], "_")[[1]]
  datalist = vector("list", length = length(snpsets))
  for (i in 1:length(snpsets)){
    variants.file = paste(variants.prefix, snpsets[i],snplabel.for.variantfile, "_matrixeqtl.traw", sep="")
    datalist[[i]] = as.data.frame(t(read.table(variants.file, sep="\t", header=T, row.names = 1)))
  }
  X = data.matrix(bind_cols(datalist))
  all.X.exclude[[j]] = which(!colnames(X.ref) %in% colnames(X))
}
cat(length(all.X.exclude), "snpsets")
for (j in seq_along(all.X.exclude)) {
  cat(all.snpsets[[j]], length(all.X.exclude[[j]]), 
      'exclude', dim(X.ref)[2] - length(all.X.exclude[[j]]), 
      "head", head(all.X.exclude[[j]]),
      "tail", tail(all.X.exclude[[j]]),'\n')
}

## Begin training
weightcol = c("gene","SNP","beta")
workingweight <- en.dir %&% "weights_bestmodel"%&% snplabel %&%".txt"
print(workingweight)
write(weightcol,file=workingweight,ncol=3,sep="\t")

set.seed(1001)

gene = colnames(expdata)[1]
for(gene in colnames(expdata)){
  cat(gene %&% '\n')
  exppheno <- expdata[,gene] ### pull expression data for gene
  #exppheno <- scale(exppheno, center=T, scale=T)  ###need to scale for fastLmPure to work properly
  #exppheno[is.na(exppheno)] <- 0
  X.exclude = all.X.exclude[[all.model.best[gene, "modelname"]]]
  lambda.min = all.model.best[gene, "lambda"]
  fit <- glmnet(X.ref,exppheno,alpha=alpha, lambda = lambda.min, exclude = X.exclude) # keep: For R2 calculation
  #plot(fit)
  saveRDS(fit, model.dir %&% "bestmodel_" %&% gene %&%".rds")
  
  ret <- as.data.frame(fit$beta[,1]) # get betas from best lambda
  ret[ret == 0.0] <- NA
  bestbetas = as.vector(ret[which(!is.na(ret)),]) # vector of non-zero betas
  names(bestbetas) = rownames(ret)[which(!is.na(ret))]
  genename <- as.character(gene)
  betafile<-cbind(genename,names(bestbetas),bestbetas) ##output "gene","SNP","beta"
  write(t(betafile),file=workingweight,ncolumns=3,append=T,sep="\t") # t() necessary for correct output from write() function
}

return()

### Generate weight files for others to use
en.dir <- "results/lasso_ukbf/"
snplabel = ".ukbflip"
workingweight <- en.dir %&% "weights_bestmodel"%&% snplabel%&%".txt"
weight = read_tsv(workingweight)
# Generate vcf file for VEP
weight$`#CHROM` = sapply(strsplit(weight$SNP,"_"), `[`, 1)
weight$POS = as.numeric(sapply(strsplit(weight$SNP,"_"), `[`, 2))
weight$REF = sapply(strsplit(weight$SNP,"_"), `[`, 3)
weight$ALT = sapply(strsplit(weight$SNP,"_"), `[`, 4)
weight %>% mutate(ID = SNP, QUAL = ".", FILTER = ".", INFO = ".") %>%
  select(`#CHROM`, POS, ID, REF, ALT, QUAL, FILTER, INFO) %>%
  distinct() %>%
  arrange(`#CHROM`, POS) %>%
  write_tsv("results/lasso_ukbf/lasso_snps_"%&% snplabel%&%".vcf")
# After VEP
weight_vep = read_tsv("results/lasso_ukbf/lasso_snps_"%&% snplabel%&%".vep.txt") %>%
  select(`#Uploaded_variation`, Existing_variation) %>% 
  rename(SNP = `#Uploaded_variation`, RSID = Existing_variation)
weight_vep = weight_vep %>% distinct()
dim(weight)
weight = weight %>% rename(weight = beta) %>% inner_join(weight_vep)
dim(weight)

# Get allele frequency
af = read_tsv("data/ukb/lasso_variant_sets/trb.hla.ukbflip.freq.afreq")
af = af %>% select(ID, ALT_FREQS) %>% rename(allele_freq = ALT_FREQS, SNP = ID)
weight.af = inner_join(weight, af)
dim(weight.af)

dir.create("models/lasso"%&% snplabel%&%"_tsv/", recursive = TRUE)
weight.af %>% 
  rename(POS_hg19 = POS) %>%
  group_split(gene) %>% 
  map(~ write_tsv(.x%>% select(`#CHROM`, POS_hg19, RSID, REF, ALT, weight, allele_freq), 
                  paste0("models/lasso"%&% snplabel%&%"_tsv/RFU", .x$gene[1]+1, "_weights.tsv")))
