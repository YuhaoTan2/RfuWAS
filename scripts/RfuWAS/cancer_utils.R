setwd("~/airr")
library(tidyverse)
library(survival)
library(broom)
library(purrr)
library(vroom)
source("scripts/ukb/pheno_utils.R")

survival_dir = "results/survival/ukbflip_hg19/"
rename = dplyr::rename
select = dplyr::select

# Diagnosis analysis ------------------------------------------
# https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=100092
# 40021	Cancer record origin 0-16
# 40019	Cancer record format 0-16
# 40009	Reported occurrences of cancer
# 40005	Date of cancer diagnosis
# 40008	Age at cancer diagnosis
# 40006	Type of cancer: ICD10
# 40013	Type of cancer: ICD9
# 40011	Histology of cancer tumour
# 40012	Behaviour of cancer tumour
# Note: https://biobank.ndph.ox.ac.uk/ukb/label.cgi?id=100074
# provides self-reported cancer
# f.34.0.0 birth year
# f.1647.0.0 birth country
# 1	England
# 2	Wales
# 3	Scotland
# 4	Northern Ireland
# 5	Republic of Ireland
# 6	Elsewhere
# -1	Do not know
# -3	Prefer not to answer


# ukb.pheno = vroom("data_xiaowei/ukb/ukb46421.tab")
# ukb.pheno.cancer = ukb.pheno %>% select(f.eid, paste0("f.40021.",0:16,".0"),
#     paste0("f.40019.",0:16,".0"), f.40009.0.0, paste0("f.40005.",0:16,".0"),
#     paste0("f.40008.",0:16,".0"), paste0("f.40006.",0:16,".0"),
#     paste0("f.40013.",0:14,".0"), paste0("f.40011.",0:16,".0"),
#     paste0("f.40012.",0:16,".0"))
# write_tsv(ukb.pheno.cancer, "data/ukb/qc/cancer.txt")
# ukb.pheno.birth = ukb.pheno %>% select(f.eid, f.34.0.0, f.1647.0.0)
# write_tsv(ukb.pheno.birth, "data/ukb/qc/birth.txt")

ukb.pheno.cancer = read_tsv("data/ukb/qc/cancer.txt")
ukb.pheno.cancer = ukb.pheno.cancer %>% rename(IID = f.eid) %>% 
  select(IID, paste0("f.40008.",0:16,".0"), paste0("f.40006.",0:16,".0"),
    paste0("f.40013.",0:14,".0"))
ukb.pheno.birth = read_tsv("data/ukb/qc/birth.txt")
ukb.pheno.birth = ukb.pheno.birth %>% rename(IID = f.eid)

#unique(ukb.pheno.cancer$f.40021.0.0)
#table((ukb.pheno.cancer %>%
#         left_join(ukb.pheno.birth) %>% filter(f.1647.0.0 == 5))$f.40021.0.0)



icd10_phecode = read_tsv("data/ukb/phenome/UKB_PHENOME_ICD10_PHECODE_MAP_20230403.txt")
icd9_phecode = read_tsv("data/ukb/phenome/UKB_PHENOME_ICD9_PHECODE_MAP_20230403.txt")

get_diagnosis_age = function(ori_phenos, phecode){
  icd10_sel = icd10_phecode %>% filter(phecode == !!phecode) %>% pull(ICD10)
  icd9_sel = icd9_phecode %>% filter(phecode == !!phecode) %>% pull(ICD9)
  icd10_sel <- gsub("\\.", "", icd10_sel)
  icd9_sel <- gsub("\\.", "", icd9_sel)

  # Gather columns into key-value pairs
  df_long <- ori_phenos %>% select(IID, paste0("f.40006.", 0:16, ".0"), paste0("f.40008.", 0:16, ".0"), paste0("f.40013.", 0:14, ".0")) %>%
    gather(key, value, -IID) %>% na.omit()

  # Separate the key column into individual columns
  df_long <- df_long %>%
    separate(key, into = c("prefix", "field", "instance", "index"), sep = "\\.")

  # Pivot the data frame to create separate columns for 40006 and 40008
  df_long <- df_long %>%
    pivot_wider(names_from = field, values_from = value) 

  # Check if columns exist, and add them if they don't
  if (!("40006" %in% colnames(df_long))) {
    df_long <- df_long %>%
      add_column("40006" = NA)
  }

  if (!("40008" %in% colnames(df_long))) {
    df_long <- df_long %>%
      add_column("40008" = NA)
  }

  if (!("40013" %in% colnames(df_long))) {
    df_long <- df_long %>%
      add_column("40013" = NA)
  }

  df_long = df_long %>%
    select(IID, `40006`, `40008`, `40013`)

  df_long_sel = df_long %>% filter(`40006` %in% icd10_sel | `40013` %in% icd9_sel)

  #df_long_sel

  # Check for inconsistent ages
  # inconsistent_age <- df_long_sel %>%
  #   group_by(IID) %>%
  #   filter(n_distinct(`40008`) > 1) %>%
  #   ungroup()

  # Print warnings for inconsistent ages
  # if (nrow(inconsistent_age) > 0) {
  #   warning("Inconsistent ages detected in the following rows: ", 
  #           paste(inconsistent_age$IID, collapse = ", "))
  # }

  #df_long_sel %>% filter(IID == "1002332")

  # Compute minimum age for each id
  result <- df_long_sel %>%
    group_by(IID) %>%
    summarise(diagnose_age = min(`40008`, na.rm = TRUE)) %>%
    ungroup()
  result
}

# # Cox regression
# # exp(coef) is the hazard ratio, coef < 0: Reduction in hazard (protective)
# # Create a function to fit the model and return the coefficients
# fit_model <- function(var_name, phenos_age_rfu) {
#   covar = c("isMale", paste0("PC", 1:10))
#   formula <- as.formula(paste0("Surv(diagnose_age, pheno) ~ ", var_name, " + ", paste(covar, collapse="+")))
#   fit <- coxph(formula, data = phenos_age_rfu)
#   tidy(fit) %>% head(1)
# }
# 
# 
prepare_phenos_age_rfu = function(phecode){
  ori_phenos = vroom(paste0("data/ukb/phenotype/X",phecode,".tsv"), show_col_types = F)
  ori_phenos = ori_phenos %>% left_join(ukb.pheno.cancer, by = "IID")
  cat("original phenos dim:", dim(ori_phenos), '\n')

  #sum(is.na(ori_phenos$f.40008.4.0))

  result = get_diagnosis_age(ori_phenos, phecode)

  # Join age column back to original data frame
  phenos_age <- ori_phenos %>% select(IID, !!paste0("X", phecode)) %>%
    left_join(result, by = "IID")
  cat("phenos dim after joining age:", dim(phenos_age), '\n')

  # Require each phenotype 1 has a corresponding age
  phenos_age = phenos_age %>%
    filter((is.na(diagnose_age) & !!sym(paste0("X", phecode)) == 0) | (!is.na(diagnose_age) & !!sym(paste0("X", phecode)) == 1))
  #df_long %>% filter(IID == "1000565")
  cat("phenos dim after requiring phenotype 1 has corresponding age:", dim(phenos_age), '\n')

  # ukb.pheno.qc = read_tsv("data/ukb/qc/sample_qc.txt")
  # ukb.pheno.qc = ukb.pheno.qc %>% rename(IID = eid) %>% select(IID, age)
  # phenos_age = phenos_age %>% left_join(ukb.pheno.qc, by = "IID")

  phenos_age = phenos_age %>% left_join(ukb.pheno.birth, by = "IID")
  # Sensoring date in
  # https://biobank.ndph.ox.ac.uk/ukb/exinfo.cgi?src=Data_providers_and_dates
  #table(phenos_age$f.1647.0.0)
  phenos_age = phenos_age %>%
    mutate(sensoring_year = 2020)
  #  mutate(sensoring_year = switch(f.1647.0.0, 1 = 2020, 2 = 2016, 3 = 2021, default = 2020))
  phenos_age = phenos_age %>% mutate(diagnose_age = ifelse(is.na(diagnose_age), sensoring_year - f.34.0.0, diagnose_age)) %>% mutate(diagnose_age = as.numeric(diagnose_age))
  cat("phenos dim after adding age for phenotype 0:", dim(phenos_age), '\n')

  sel_expression = ori_expression %>% filter(IID %in% phenos_age$IID)

  phenos_age_rfu = phenos_age %>% rename(pheno = !!paste0("X",phecode)) %>%
    left_join(sel_expression, by = "IID")

  phenos_age_rfu = phenos_age_rfu %>% left_join(sample_qc, by = "IID")

  cat("phenos dim after adding rfu and qc data:", dim(phenos_age_rfu), '\n')
  phenos_age_rfu
}
# 
# 
# phecode_cox = function(phecode){
#   cat <- function(...) base::cat(..., file = "log/cancer.out", append = TRUE)
#   cat(phecode, '\n')
#   # get_description(phecode)
#  
#   phenos_age_rfu = prepare_phenos_age_rfu(phecode)
# 
#   # Apply the function to each column name and row bind the results
#   cox_coef <- map2_dfr(all_rfu, list(phenos_age_rfu), fit_model)
#   cox_coef = cox_coef %>% arrange(p.value)
#   write_tsv(cox_coef, paste0("results/coxph/", phecode, ".txt"))
# }

# Survival analysis -------------------------------
# 40023	Records in death dataset
# 40018	Death record format
# 40020	Death record origin
# 40000	Date of death
# 40007	Age at death
# 40001	Underlying (primary) cause of death: ICD10
# 40002	Contributory (secondary) causes of death: ICD10
# 40010	Description of cause of death

# ukb.pheno = vroom("data_xiaowei/ukb/ukb46421.tab")
# ukb.pheno.death = ukb.pheno %>% select(f.eid, paste0("f.40018.",0:1,".0"),
#     paste0("f.40020.",0:1,".0"), paste0("f.40000.",0:1,".0"),
#     paste0("f.40007.",0:1,".0"), paste0("f.40001.",0:1,".0"),
#     paste0("f.40002.0.",1:14), paste0("f.40002.1.",1:14),
#     "f.40010.0.0")
# write_tsv(ukb.pheno.death, "data/ukb/qc/death.txt")

ukb.pheno.death = read_tsv("data/ukb/qc/death.txt")
ukb.pheno.death = ukb.pheno.death %>% rename(IID = f.eid) %>% 
  select(IID, paste0("f.40007.",0:1,".0"), paste0("f.40001.",0:1,".0")) 
ukb.pheno.death = ukb.pheno.death %>% 
  mutate(death_age = pmin(f.40007.0.0, f.40007.1.0, na.rm = T), 
         death = ifelse(is.na(death_age), 0, 1))

# whole_phenos = vroom("data/ukb/phenome/UKB_PHENOME_NO_EXCLUSIONS_20230403.txt")
# cancer_phenos = whole_phenos %>% select(IID, all_of(intersect(paste0("X", cancer_phecode), colnames(.))))
# cancer_phenos %>% write_tsv("data/ukb/phenome/cancer_phenos_no_exclusions.txt")

cancer_phenos = read_tsv("data/ukb/phenome/cancer_phenos_no_exclusions.txt")

prepare_phenos_death_rfu = function(phecode, permute = F){
  ori_phenos = cancer_phenos %>% select(IID, !!sym(paste0("X", phecode))) 
  ori_phenos_sel = ori_phenos %>% filter(!!sym(paste0("X", phecode)) == 1)
  phenos_death = ukb.pheno.death %>% 
    filter(IID %in% ori_phenos_sel$IID) %>%
    left_join(ukb.pheno.cancer, by = "IID")
  cat("death infomation with disease dim:", dim(phenos_death), '\n')

  result = get_diagnosis_age(phenos_death, phecode)
  
  # Join age column back to original data frame
  # End-point: all-cause mortality
  phenos_death <- phenos_death %>% 
    select(IID, death_age, death) %>% inner_join(result, by = "IID")
  cat("phenos dim after joining age:", dim(phenos_death), '\n')

  phenos_death = phenos_death %>% left_join(ukb.pheno.birth, by = "IID")
  # Sensoring date in
  # https://biobank.ndph.ox.ac.uk/ukb/exinfo.cgi?src=Data_providers_and_dates
  phenos_death = phenos_death %>% 
    mutate(sensoring_year = 2022) %>% 
    mutate(death_age = ifelse(is.na(death_age), sensoring_year - f.34.0.0, death_age), diagnose_age = as.numeric(diagnose_age)) %>% 
    mutate(survive_year = death_age - diagnose_age)
  if (permute){
    print(phenos_death)
    print("Permuting ...")
    phenos_death <- phenos_death %>% 
               mutate(IID = sample(IID))
    print(phenos_death)
  }
  cat("phenos dim after adding age for phenotype 0:", dim(phenos_death), '\n')
  
  phenos_death_rfu = phenos_death %>% 
    inner_join(ori_expression, by = "IID") %>% inner_join(sample_qc, by = "IID")
  
  cat("phenos dim after adding rfu and qc data:", dim(phenos_death_rfu), '\n')
  
  return(phenos_death_rfu)
}

# Cox regression
# exp(coef) is the hazard ratio, coef < 0: Reduction in hazard (protective)
# Create a function to fit the model and return the coefficients
fit_survival_model <- function(var_name, phenos_death_rfu, covar) {
  phenos_death_rfu_sel = phenos_death_rfu %>% select(survive_year, death, all_of(var_name), all_of(covar)) %>% na.omit()
  true_num = phenos_death_rfu_sel %>% filter(death == 1) %>% nrow()
  false_num = phenos_death_rfu_sel %>% filter(death == 0) %>% nrow()
  formula <- as.formula(paste0("Surv(survive_year, death) ~ ", var_name, " + ", paste(covar, collapse="+")))
  if (true_num == 0 | false_num == 0) return(data.frame())
  tryCatch({
    fit <- coxph(formula, data = phenos_death_rfu_sel, iter.max = 1000)
    return(tidy(fit) %>% head(1) %>% mutate(true_num = !!true_num, false_num = !!false_num))
  }, error = function(e) {
    cat("**Error in ", phenos_death_rfu %>% head(1) %>% pull(phecode), var_name, conditionMessage(e), '\n')
    return(data.frame())
  }, warning = function(w) {
    cat("**Warning in ", phenos_death_rfu %>% head(1) %>% pull(phecode), var_name, conditionMessage(w), '\n')
    fit <- coxph(formula, data = phenos_death_rfu_sel)
    return(tidy(fit) %>% head(1) %>% mutate(true_num = !!true_num, false_num = !!false_num))
  })
}

phecode_survival = function(phecode, permute = F){
  log_file = paste0("log/cancer_survival.out")
  #cat <- function(...) base::cat(..., file = log_file, append = TRUE)
  #sink(log_file, append = TRUE)
  cat(phecode, '\n')
  # get_description(phecode)
 
  phenos_death_rfu = prepare_phenos_death_rfu(phecode, permute = permute) %>% mutate(phecode = !!phecode)
  true_num = phenos_death_rfu %>% filter(death == 1) %>% nrow()
  false_num = phenos_death_rfu %>% filter(death == 0) %>% nrow()

  # As suggested in https://www.nature.com/articles/6601120
  # at least 10 events need to be observed for each covariate considered
  # We have 13 variates, so we require 130 events.

  if (true_num == 0 | false_num == 0){
    cat(phecode, dim(phenos_death_rfu), " No death or no survival, skip\n")
    return()
  }

  # Apply the function to each column name and row bind the results
  if (true_num >= 130){
    cat(phecode, true_num, "Calculating 13 cov\n")
    cox_coef <- pmap_dfr(list(
      var_name = all_rfu, phenos_death_rfu = list(phenos_death_rfu), covar = list(c("diagnose_age", "isMale", paste0("PC", 1:10)))), 
      fit_survival_model)
    cox_coef = cox_coef %>% arrange(p.value)
    if (permute) {
      write_tsv(cox_coef, paste0(survival_dir, phecode, "_survival_13cov_permute.txt"))
    } else {
      write_tsv(cox_coef, paste0(survival_dir, phecode, "_survival_13cov.txt"))
    }
  }
}
