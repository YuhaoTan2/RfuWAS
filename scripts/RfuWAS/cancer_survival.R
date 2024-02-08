setwd("~/airr")
library(tidyverse)
library(survival)
library(broom)
library(purrr)
library(vroom)
source("scripts/RfuWAS/cancer_utils.R")
source("scripts/RfuWAS/ukb_utils.R")


if(!dir.exists(survival_dir)) dir.create(survival_dir)


set.seed(1234)
for(phecode in cancer_phecode) phecode_survival(phecode, permute = T)
