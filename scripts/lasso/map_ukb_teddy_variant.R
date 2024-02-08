### Find the intersection between UKB and TEDDY variants after mapping to the same reference genome
### If the variant name is different, swap the ref and alt allele
### Output the variant name in TEDDY and UKB, and the mapping between them

setwd("~/airr")
"%&%" = function(a,b) paste(a,b,sep="")
library(tidyverse)
library(readr)
library(purrr)
source("scripts/ukb/ukb_utils.R")

# UKB -----------------------------------------------------------------
teddy.assembly = "hg19"
if (teddy.assembly == "hg38"){
  teddy.var = read_tsv("data/gwas/sel_trb_hla.bim", col_names = FALSE)$V1
  output_prefix = "data/ukb/genotype/"
} else {
  teddy.var = read_tsv("data/gwas/sel_trb_hla_hg19.pvar")$ID
  output_prefix = "data/ukb/qc/"
}
cat("TEDDY variants number:", length(teddy.var), '\n')

ukb.impute = T
if (ukb.impute){
  ukb.pvar = read_tsv("data/ukb/qc/ukb_imp_hwe.pvar")$ID
} else {
  ukb.pvar = read_tsv("data/ukb/genotype/sel_trb_hla.pvar")$ID
}

cat("UKB variants number:", length(ukb.pvar), '\n')

## Find variants exist in both UKB and TEDDY
# More thorough way
matches = map_df1_df2(teddy.var, ukb.pvar)

# UKB name
write.table(matches$Name_df2, paste0(output_prefix, "ukb.varid"), sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote=F)
write.table(matches$Name_df1, paste0(output_prefix, "teddy.varid"), sep = "\t",
            col.names = FALSE, row.names = FALSE, quote=F)
# Convert UKB name to TEDDY name
write.table(matches %>% select(Name_df2, Name_df1), paste0(output_prefix, "ukb.old.new.varid"), sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote=F)
# ref allele file
write.table(matches %>% select(Name_df1, REF), paste0(output_prefix, "ukb.ref"), sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote=F)

# We know that sometimes swapping ref and alt may align different variants, so we validate the intersection by changing the reference
teddy.pvar = read_tsv("~/airr/data/gwas/sel_trb_hla_hg19.pvar")
teddy.pvar$POS = sapply(strsplit(teddy.pvar$ID, "_"), function(x) x[2])
write_tsv(teddy.pvar, "~/airr/data/gwas/sel_trb_hla_hg19_repos.pvar")
file.copy("~/airr/data/gwas/sel_trb_hla_hg19.pgen", "~/airr/data/gwas/sel_trb_hla_hg19_repos.pgen")
file.copy("~/airr/data/gwas/sel_trb_hla_hg19.psam", "~/airr/data/gwas/sel_trb_hla_hg19_repos.psam")

ukb.pvar.ref = read_tsv("~/airr/data/ukb/qc/ukb_imp_hwe_ref.pvar", comment = "##")
ukb.pvar.ref = paste(ukb.pvar.ref[["#CHROM"]], ukb.pvar.ref$POS, ukb.pvar.ref$REF, ukb.pvar.ref$ALT, sep = "_")
teddy.pvar.ref = read_tsv("~/airr/data/gwas/sel_trb_hla_hg19_repos_ref.pvar", comment = "##")
teddy.pvar.ref = paste(teddy.pvar.ref[["#CHROM"]], teddy.pvar.ref$POS, teddy.pvar.ref$REF, teddy.pvar.ref$ALT, sep = "_")
teddy.ukb.ref = intersect(teddy.pvar.ref, ukb.pvar.ref)
teddy.ukb.ref.reverse = gsub("^(\\d+)_(\\d+)_(.+)_(.+)$", "\\1_\\2_\\4_\\3", teddy.ukb.ref)
setdiff(ukb.name.intersect, c(teddy.ukb.ref, teddy.ukb.ref.reverse))

# Emerson -----------------------------------------------------
teddy_ukb.pvar = read_tsv("data/ukb/lasso_variant_sets/ukb_intersect_teddyname_teddyallele.pvar")$ID
emerson.pvar = read_tsv("data/STAMPEED/trb_hla/STAMPEED.trb.hla.pvar")$ID
length(teddy_ukb.pvar)
length(emerson.pvar)

## Find variants exist in both Emerson and TEDDY
# More thorough way
matches = map_df1_df2(teddy_ukb.pvar, emerson.pvar)

# teddy_ukb name
output_prefix = "data/STAMPEED/trb_hla/"
write.table(matches$Name_df2, paste0(output_prefix, "emerson.varid"), sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote=F)
write.table(matches$Name_df1, paste0(output_prefix, "teddy_ukb.varid"), sep = "\t",
            col.names = FALSE, row.names = FALSE, quote=F)
# Convert UKB name to TEDDY name
write.table(matches %>% select(Name_df2, Name_df1), paste0(output_prefix, "emerson.old.new.varid"), sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote=F)
# ref allele file
write.table(matches %>% select(Name_df1, REF), paste0(output_prefix, "emerson.ref"), sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote=F)

##########################
# Overlap with not imputed variants
# ukb.hg38.var = read_tsv("data/ukb/genotype/ukb_intersect_teddyname_teddyallele.pvar")
# teddy.hg38.to.hg19 = read_tsv("data/gwas/teddy.old.new.varid", col_names = F)
# teddy.hg38.to.hg19 = teddy.hg38.to.hg19 %>% rename(ID = X1, ID.hg19 = X2)
# dim(ukb.hg38.var)
# ukb.hg19.var = inner_join(ukb.hg38.var, teddy.hg38.to.hg19)
# ukb.hg19.pvar = ukb.hg19.var$ID.hg19

