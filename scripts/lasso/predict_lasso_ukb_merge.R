setwd("~/airr")
"%&%" = function(a,b) paste0(a,b)
library(tidyverse)
library(glmnet)
library(vroom)

dataset.name = "ukb_hg19"
snplabel = ".ukbflip"

# rfu.df: rfus * samples
# pred_exp: samples * rfus
file_paths <- paste0("results/lasso_ukbf/", "bestmodel", snplabel, "_", dataset.name, "_", 1:24, ".txt")
rfu.df <- map(file_paths, vroom)
rfu.normalized <- bind_cols(rfu.df[[1]], !!!map(rfu.df[-1], ~ .x[, -1]))
dim(rfu.normalized)
pred_exp = as.data.frame(t(rfu.normalized))
colnames(pred_exp) = pred_exp[1,]
pred_exp = pred_exp[-1,]
dim(pred_exp)
write_tsv(pred_exp, paste0("results/lasso_ukbf/", "bestmodel",snplabel,"_", dataset.name, "_merge.txt"))
writeLines(rownames(pred_exp), paste0("results/lasso_ukbf/", "bestmodel",snplabel,"_", dataset.name, "_merge.txt", ".rownames"))
