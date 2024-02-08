# QC information of samples that pass QC is saved to data/ukb/qc/sample_qc.kt
# Sample ID that passes QC is saved to data/ukb/qc/sample_filtered.id
# Phenome after QC is saved to data/ukb/qc/UKB_PHENOME_DESCRIPTION_QC.txt

setwd("~/airr")
library(vroom)
library(tidyverse)

### 1. Sample QC
# 1.1 Filter samples based on phenotype file
ukb.pheno = vroom("data_xiaowei/ukb/ukb46421.tab")
ukb.pheno.qc = ukb.pheno %>% select(f.eid, f.21003.0.0, f.22001.0.0,
                                    f.22019.0.0, f.22020.0.0, f.22006.0.0,
                                    f.22021.0.0, f.22009.0.1, f.22009.0.2,
                                    f.22009.0.3, f.22009.0.4, f.22009.0.5,
                                    f.22009.0.6, f.22009.0.7, f.22009.0.8,
                                    f.22009.0.9, f.22009.0.10)
ukb.pheno.qc = ukb.pheno.qc %>% 
  rename(eid = f.eid, age = f.21003.0.0, Inferred.Gender = f.22001.0.0,
         putative.sex.chromosome.aneuploidy = f.22019.0.0,
         used.in.pca.calculation = f.22020.0.0,
         in.white.British.ancestry.subset = f.22006.0.0,
         excess.relatives = f.22021.0.0,
         PC1 = f.22009.0.1, PC2 = f.22009.0.2,
         PC3 = f.22009.0.3, PC4 = f.22009.0.4, PC5 = f.22009.0.5,
         PC6 = f.22009.0.6, PC7 = f.22009.0.7, PC8 = f.22009.0.8,
         PC9 = f.22009.0.9, PC10 = f.22009.0.10)

table(ukb.pheno.qc$Inferred.Gender)
sum(is.na(ukb.pheno.qc$Inferred.Gender))
table(ukb.pheno.qc$putative.sex.chromosome.aneuploidy)
table(ukb.pheno.qc$used.in.pca.calculation)
table(ukb.pheno.qc$in.white.British.ancestry.subset)
table(ukb.pheno.qc$excess.relatives)

# Before filtering
write_tsv(ukb.pheno.qc, "data/ukb/qc/sample_qc.txt")

ukb.pheno.qc %>% filter(eid < 0)
ukb.filter = ukb.pheno.qc %>% filter(in.white.British.ancestry.subset == 1)
ukb.filter = ukb.filter %>% filter(used.in.pca.calculation == 1)
ukb.filter = ukb.filter %>%
  filter(is.na(putative.sex.chromosome.aneuploidy))

# 1.2 Filter samples with genotype data
fam = vroom("data_from_xiaowei/ukb/bfiles_all_imputed/ukb_imp_chr6_v3.6.27477797.34448354.fam",
            col_names = NULL)
sampleid = fam$X1
ukb.filter = ukb.filter %>% filter(eid %in% sampleid)

ukb.filter = ukb.filter %>% 
  mutate(sample = eid, isMale = Inferred.Gender) %>%
  select(sample, isMale, PC1, PC2, PC3, 
         PC4, PC5, PC6, PC7, PC8, PC9, PC10, age)
# After filtering
write_tsv(ukb.filter, "data/ukb/qc/sample_qc.kt")
ukb.filter.sample = ukb.filter %>% 
  mutate(`#FID` = sample, IID = sample) %>%
  select(`#FID`, IID)
write_tsv(ukb.filter.sample, "data/ukb/qc/sample_filtered.id")

### 2. Genotype QC
# 2.1 Variants with INFO > 0.8
all.variants = vroom("data/ukb/qc/ukb_imp_trb_hla.pvar")
all.variants.id = all.variants$ID

info.chr6 = vroom("data/ukb/qc/ukb_mfi_chr6_v3.txt", col_names = NULL)
names(info.chr6) = c("Alternate_id", "RS_id", "Position", 
                     "Allele1", "Allele2", "MAF", "Minor Allele",
                     "Info score")
info.chr6 = info.chr6 %>% filter(Position >= 28477797-1e6 & Position <= 33448354+1e6)
info.chr6 = info.chr6 %>% mutate(chr = 6)

info.chr7 = vroom("data/ukb/qc/ukb_mfi_chr7_v3.txt", col_names = NULL)
names(info.chr7) = c("Alternate_id", "RS_id", "Position", 
                     "Allele1", "Allele2", "MAF", "Minor Allele",
                     "Info score")
info.chr7 = info.chr7 %>% filter(Position >= 141998851-1e6 & Position <= 142510972+1e6)
info.chr7 = info.chr7 %>% mutate(chr = 7)

info.chr67 = bind_rows(info.chr6, info.chr7)
info.chr67 = info.chr67 %>%
  mutate(forward.ID = paste(chr, Position, Allele1, Allele2, sep='_'),
         reverse.ID = paste(chr, Position, Allele2, Allele1, sep='_'))

info.filter = info.chr67 %>% filter(`Info score` > 0.8)
dim(info.chr67)
dim(info.filter)
info.filter.ID = c(info.filter$forward.ID, info.filter$reverse.ID)
info.filter.ID = intersect(info.filter.ID, all.variants.id)
length(info.filter.ID)
writeLines(info.filter.ID, "data/ukb/qc/variants.pass.info.rsid")

### 3. Phenotype QC, generate phenotype file for each phenotype
ukb.pheno = vroom("data/ukb/phenome/UKB_PHENOME_20230403.txt")
ukb.filter = read_tsv("data/ukb/qc/sample_qc.kt")
ukb.pheno.f = ukb.pheno %>% filter(IID %in% ukb.filter$sample)
ukb.pheno.f = as.data.frame(ukb.pheno.f)
rownames(ukb.pheno.f) = ukb.pheno.f$IID
ukb.pheno.f = ukb.pheno.f[, -1]

# calculate the ratio of 1 over (0 or 1) for each column
n_ones <- colSums(ukb.pheno.f == 1, na.rm = TRUE)
n_zeros_ones <- colSums(!is.na(ukb.pheno.f), na.rm = TRUE)
ratios <- n_ones / n_zeros_ones
cat(sum(ratios == 0), "phenotypes have no non-missing values\n")

# subset the data frame based on the ratio threshold
ukb.pheno.filtered <- ukb.pheno.f[, ratios > 0.001]
ukb.pheno.notpass = ukb.pheno.f[, ratios <= 0.001]
ukb.pheno.des = vroom("data/ukb/phenome/UKB_PHENOME_DESCRIPTION_20230403.txt")

save_pheno_qc = function(ukb.pheno.filtered, filter_name){
  dim(ukb.pheno.filtered)
  ukb.pheno.filtered$IID = rownames(ukb.pheno.filtered)
  write_tsv(ukb.pheno.filtered, paste0("data/ukb/qc/UKB_PHENOME_",filter_name,".txt"))

  ukb.pheno.des.f = ukb.pheno.des %>% 
    filter(phecode %in% gsub("X", "", colnames(ukb.pheno.filtered)))
  write_tsv(ukb.pheno.des.f, paste0("data/ukb/qc/UKB_PHENOME_DESCRIPTION_",filter_name,".txt"))

  if (!file.exists("data/ukb/phenotype/")) {
    dir.create("data/ukb/phenotype/")
  }
  for (col in names(ukb.pheno.filtered)[-length(names(ukb.pheno.filtered))]) {
    col_data = ukb.pheno.filtered %>% select(IID, all_of(col))
    write_tsv(col_data, file = paste0("data/ukb/phenotype/", col, ".tsv"))
  }
}

save_pheno_qc(ukb.pheno.filtered, "QC")
save_pheno_qc(ukb.pheno.notpass, "NOTPASS")

set.seed(1234)
ukb.pheno.filtered = read_tsv(paste0("data/ukb/qc/UKB_PHENOME_","QC",".txt"))
for (col in names(ukb.pheno.filtered)[-length(names(ukb.pheno.filtered))]) {
  col_data = ukb.pheno.filtered %>% select(IID, all_of(col)) %>% mutate(IID = sample(IID))
  write_tsv(col_data, file = paste0("data/ukb/phenotype/", col, "_permuted.tsv"))
}

