library(tidyverse)
library(vroom)

stratify_data_fn = "data/ukbflip/stratify_tmp"
stratify_result_fn = "results/assoc_pheno/ukbflip_hla_stratify"

# Get all neoplasm phecode
phenome_des = vroom("data/ukb/qc/UKB_PHENOME_DESCRIPTION_QC.txt", show_col_types = FALSE)
phenome_des_notpass = vroom("data/ukb/qc/UKB_PHENOME_DESCRIPTION_NOTPASS.txt", show_col_types = FALSE)
phenome_des_all = rbind(phenome_des, phenome_des_notpass)
cancer_phecode = phenome_des_all %>% filter(group == "neoplasms") %>% pull(phecode)

get_rfu_num = function(){
  output_dir = "results/assoc_pheno/ukbflip_hg19/"
  assoc.fn = paste(output_dir,"assoc_X", "174.1",".tsv", sep = "")
  assoc_df = read_tsv(assoc.fn, show_col_types = FALSE)
  rfu.num = nrow(assoc_df)
  rfu.num
}

get_predictable_rfu = function(){
  output_dir = "results/assoc_pheno/ukbflip_hg19/"
  assoc.fn = paste(output_dir,"assoc_X", "174.1",".tsv", sep = "")
  assoc_df = read_tsv(assoc.fn, show_col_types = FALSE)
  assoc_df$gene
}

get_description = function(phecode_input){
    phenome_des_all %>% filter(phecode == phecode_input) %>% pull(description)
}

get_phecode = function(des_input){
    phenome_des_all %>% filter(str_detect(des_input, description)) %>% select(phecode, description)
}


map_df1_df2 = function(df1.var, df2.var){
  df1 <- tibble(Name_df1 = unlist(df1.var)) %>%
    separate(Name_df1, into = c("CHR", "POS", "REF", "ALT"), sep = "_", remove = FALSE)
  df2 <- tibble(Name_df2 = unlist(df2.var)) %>%
    separate(Name_df2, into = c("CHR", "POS", "REF", "ALT"), sep = "_", remove = FALSE)

  # Function to flip nucleotides
  flip_nucleotide <- function(sequence) {
    flipped_sequence <- chartr("ATCG", "TAGC", sequence)
    return(flipped_sequence)
  }

  # Apply the function to flip nucleotides in df1
  df1_flipped <- df1 %>%
    mutate(across(c(REF, ALT), ~map_chr(., flip_nucleotide)))

  # Initialize an empty tibble to store matches
  matches <- tibble(Name_df1 = character(), Name_df2 = character(), REF = character(), type = character())

  # Match variants with REF, ALT matched
  matched <- inner_join(df1, df2, by=c("CHR", "POS", "REF", "ALT"))
  if (nrow(matched) > 0) { # 24133
    matches <- bind_rows(matches, tibble(Name_df1 = matched$Name_df1, Name_df2 = matched$Name_df2, REF = matched$REF, type = "matched"))
    df1 <- filter(df1, !Name_df1 %in% matched$Name_df1)
    df2 <- filter(df2, !Name_df2 %in% matched$Name_df2)
  }

  # Match with REF ALT reversed
  reversed <- inner_join(df1, df2, by=c("CHR" = "CHR", "POS" = "POS", "REF" = "ALT", "ALT" = "REF"))
  if (nrow(reversed) > 0) { # 6001
    matches <- bind_rows(matches, tibble(Name_df1 = reversed$Name_df1, Name_df2 = reversed$Name_df2, REF = reversed$REF, type = "reversed"))
    df1 <- filter(df1, !Name_df1 %in% reversed$Name_df1)
    df2 <- filter(df2, !Name_df2 %in% reversed$Name_df2)
  }

  # Match with REF ALT flipped
  flipped <- inner_join(df1_flipped, df2, by=c("CHR", "POS", "REF", "ALT"))
  if (nrow(flipped) > 0) { # 159
    matches <- bind_rows(matches, tibble(Name_df1 = flipped$Name_df1, Name_df2 = flipped$Name_df2, REF = flipped$REF, type = "flipped"))
    df1 <- filter(df1, !Name_df1 %in% flipped$Name_df1)
    df2 <- filter(df2, !Name_df2 %in% flipped$Name_df2)
  }

  # Match with REF ALT flipped and reversed
  flipped_reversed <- inner_join(df1_flipped, df2, by=c("CHR" = "CHR", "POS" = "POS", "REF" = "ALT", "ALT" = "REF"))
  if (nrow(flipped_reversed) > 0) { # 72
    matches <- bind_rows(matches, tibble(Name_df1 = flipped_reversed$Name_df1, Name_df2 = flipped_reversed$Name_df2, REF = flipped_reversed$REF, type = "flipped_reversed"))
    df1 <- filter(df1, !Name_df1 %in% flipped_reversed$Name_df1)
    df2 <- filter(df2, !Name_df2 %in% flipped_reversed$Name_df2)
  }

  print(dim(matched))
  print(dim(reversed))
  print(dim(flipped))
  print(dim(flipped_reversed))
  return(matches)
}

# calculating association after HLA stratification --------------------
extract_hla = function(){
    # Extract HLA allele information
    ukb.pheno = vroom("data_xiaowei/ukb/ukb46421.tab")
    ukb.pheno.hla = ukb.pheno %>% select(f.eid, f.22182.0.0)
    dim(ukb.pheno.hla)
    ukb.pheno.hla = ukb.pheno.hla %>% filter(!is.na(f.22182.0.0))
    dim(ukb.pheno.hla) # After filtering NA: 488239
    ukb.pheno.hla$comb = paste(ukb.pheno.hla$f.eid, ukb.pheno.hla$f.22182.0.0, sep = ",")
    write_tsv(ukb.pheno.hla %>% select(comb), "data/ukb/qc/hla.txt", col_names = FALSE)
}


prepare_samples = function(sample, sample_name, phecode, ori_phenos = ori_phenos, printcommand = T){
    # sample_name: file name
    # Find the position of sample in the original phenotype file
    # if sample is data.frame
    if (is.data.frame(sample)){
        sample = sample$f.eid
    }
    idx = match(sample, ori_phenos$IID)
    idx = idx[!is.na(idx)]
    
    phenos = ori_phenos[idx,]
    phenos.table = table(phenos %>% pull(paste0("X", phecode)))
    print(phenos.table)
    if ( ( (phenos.table["0"] < 5) | (is.na(phenos.table["0"])) ) | 
      ( (phenos.table["1"] < 5) | (is.na(phenos.table["1"])) ) ) {
      print(paste0(phecode, "has < 5 true or negative."))
      return(F)}
    
    expression = ori_expression[idx,]
    covariats = ori_covariates[idx,]

    # if (file.exists(paste0("",stratify_data_fn,"/X",phecode,"_", sample_name, ".tsv"))) return()

    write_tsv(phenos, paste0("",stratify_data_fn,"/X",phecode,"_", sample_name, ".tsv"))
    if (!file.exists(paste0("",stratify_data_fn,"/expr_",sample_name,".txt"))){
        write_tsv(expression, paste0("",stratify_data_fn,"/expr_",sample_name,".txt"))
    }
    if (!file.exists(paste0("",stratify_data_fn,"/cov_",sample_name,".txt"))){
        write_tsv(covariats, paste0("",stratify_data_fn,"/cov_",sample_name,".txt"))
    }
    
    if (printcommand) print(paste0("python3 ~/airr/othercode/MetaXcan/software/PrediXcanAssociation.py ",
        "--expression_file ~/airr/",stratify_data_fn,"/expr_",sample_name,".txt ",
        "--input_phenos_file ~/airr/",stratify_data_fn,"/X",phecode,"_", sample_name, ".tsv ",
        "--input_phenos_column X", phecode, " ",
        "--covariates_file ~/airr/",stratify_data_fn,"/cov_",sample_name,".txt ",
        "--covariates isMale PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 age ",
        "--output ~/airr/",stratify_result_fn,"/assoc_", phecode, "_", sample_name, ".tsv ", 
        "--mode logistic"))
    return(T)
}

# Function to filter and prepare samples
filter_and_prepare_samples <- function(ukb_pheno_hla, condition, name, phecode) {
  select = dplyr::select
  # Filter based on HLA condition
  # Filter samples
  sample = ukb_pheno_hla %>%
    filter(eval(rlang::parse_expr(condition))) %>%
    select(f.eid)
  
  # Not sample
  not_sample = ukb_pheno_hla %>%
    filter(!eval(rlang::parse_expr(condition))) %>%
    select(f.eid)
  
  # Assert that the counts add up to the total number of samples
  assertthat::assert_that(length(sample$f.eid) + length(not_sample$f.eid) == length(ukb_pheno_hla$f.eid))
  
  # Prepare for PrediXcan input
  prepare_samples(sample, name, phecode) 
  prepare_samples(not_sample, paste0("not",name), phecode)
}

find_sig_validate = function(phecode, allele1, allele2, validate.sig, rfu_list = 0:4999){
  select = dplyr::select
  assoc.full = vroom(paste0("results/assoc_pheno/ukbflip_hg19/assoc_X", phecode, ".tsv"))
  assoc.allele1 = vroom(paste0("",stratify_result_fn,"/assoc_", phecode, "_", allele1, ".tsv"))
  assoc.allele2 = vroom(paste0("",stratify_result_fn,"/assoc_", phecode, "_", allele2, ".tsv"))
  assoc.notallele1 = vroom(paste0("",stratify_result_fn,"/assoc_", phecode, "_not", allele1, ".tsv"))
  assoc.notallele2 = vroom(paste0("",stratify_result_fn,"/assoc_", phecode, "_not", allele2, ".tsv"))

  pvalue.thre = 0.05 / nrow(assoc.full)
  full.gene = assoc.full %>% filter(pvalue < pvalue.thre, gene %in% rfu_list) %>% select(gene)
  print("full")
  print(dim(full.gene))

  allele1.gene = assoc.allele1 %>% filter(pvalue < pvalue.thre, gene %in% rfu_list) %>% select(gene)
  print(allele1)
  print(dim(allele1.gene))
  notallele1.gene = assoc.notallele1 %>% filter(pvalue < pvalue.thre, gene %in% rfu_list) %>% select(gene)
  print(dim(notallele1.gene))

  allele2.gene = assoc.allele2 %>% filter(pvalue < pvalue.thre, gene %in% rfu_list) %>% select(gene)
  print(allele2)
  print(dim(allele2.gene))
  notallele2.gene = assoc.notallele2 %>% filter(pvalue < pvalue.thre, gene %in% rfu_list) %>% select(gene)
  print(dim(notallele2.gene))

  full.u.pvalue = validate.sig %>% filter(RFU %in% full.gene$gene)
  print("full")
  #print(full.u.pvalue)
  print(length(unique(full.u.pvalue$RFU)))

  allele1.u.pvalue = validate.sig %>% filter(RFU %in% allele1.gene$gene)
  print(allele1)
  #print(allele1.u.pvalue)
  print(length(unique(allele1.u.pvalue$RFU)))
  notallele1.u.pvalue = validate.sig %>% filter(RFU %in% notallele1.gene$gene)
  #print(notallele1.u.pvalue)
  print(length(unique(notallele1.u.pvalue$RFU)))

  allele2.u.pvalue = validate.sig %>% filter(RFU %in% allele2.gene$gene)
  print(allele2)
  #print(allele2.u.pvalue)
  print(length(unique(allele2.u.pvalue$RFU)))
  notallele2.u.pvalue = validate.sig %>% filter(RFU %in% notallele2.gene$gene)
  #print(notallele2.u.pvalue)
  print(length(unique(notallele2.u.pvalue$RFU)))
  # return: 1 - 5: overlap RFU. 6 - 10: all significant genes
  return (list(full.u.pvalue, allele1.u.pvalue, notallele1.u.pvalue, allele2.u.pvalue, notallele2.u.pvalue, full.gene, allele1.gene, notallele1.gene, allele2.gene, notallele2.gene))
}


sig_num_pvalue = function(all.gene, sel.rfu.num, overlap.num, validate.sig.rfu){
  intersect.nums = c()
  total_sample_num = 10000
  for (i in 1:total_sample_num){
    sel.rfu = sample(all.gene, sel.rfu.num, replace = FALSE)
    intersect.num = length(intersect(sel.rfu, validate.sig.rfu))
    intersect.nums = c(intersect.nums, intersect.num)
  }
  p = ggplot(data = data.frame(intersect.nums = intersect.nums), aes(x = intersect.nums)) +
    geom_histogram(fill = "steelblue", color = "white", binwidth = 1) +
    geom_vline(xintercept = overlap.num, color = "green", linetype = "dashed") +
    labs(title = "Histogram", x = "Validated number", y = "Frequency")
  print(p)
  print(sum(intersect.nums >= overlap.num) / total_sample_num)
  return(intersect.nums)
}

get_sig = function(phecode, allele1){
  assoc.allele1 = vroom(paste0("",stratify_result_fn,"/assoc_", phecode, "_", allele1, ".tsv"))
  pvalue.thre = 0.05 / nrow(assoc.allele1)
  add_dq2.5.gene = assoc.allele1 %>% filter(pvalue < pvalue.thre) %>% dplyr::select(gene) %>% pull()
  print(length(add_dq2.5.gene))
  return(add_dq2.5.gene)
}

plot.venn = function(set1, set2, set3, set1.name, set2.name, set3.name){
  sets <- list(set1 = set1, set2 = set2, set3 = set3)
  temp = venn.diagram(x = sets, category.names = c(set1.name, set2.name, set3.name) ,filename = NULL)
  grid.newpage()
  grid.draw(temp)
}

# calculating association after disease/variants stratification --------------------
pathogenic_gene_score = function(ukb_gene_var_fn, score_coef_fn){
  # Create score coef file for score calculation
  ukb_NF1_var = read_tsv(ukb_gene_var_fn, comment = "##")
  ukb_NF1_vid = ukb_NF1_var$ID

  score_coef = data.frame(
    variant = ukb_NF1_vid,
    allele = sapply(strsplit(ukb_NF1_vid, "_"), function(x) x[[4]]),
    score = rep(1, length(ukb_NF1_vid))
  )
  print(head(score_coef))
  print(dim(score_coef))
  write_tsv(score_coef, score_coef_fn)
}

cal_true_false_for_disease = function(disease_IID){
  disease_IID_num = length(disease_IID)
  all_phecodes = c(phenome_des$phecode, phenome_des_notpass$phecode)
  all_phenos_true_num = c()
  all_phenos_false_num = c()
  for (current_phecode in all_phecodes) {
    ori_phenos = vroom(paste0("data/ukb/phenotype/X",current_phecode,".tsv"), show_col_types = FALSE)
    idx = match(disease_IID, ori_phenos$IID)
    idx = idx[!is.na(idx)]
    phenos = ori_phenos[idx,]
    phenos_true_num = phenos %>% filter(!!sym(paste0("X",current_phecode)) == 1) %>% pull()
    phenos_true_num = length(phenos_true_num)
    phenos_false_num = phenos %>% filter(!!sym(paste0("X",current_phecode)) == 0) %>% pull()
    phenos_false_num = length(phenos_false_num)
    all_phenos_true_num = c(all_phenos_true_num, phenos_true_num)
    all_phenos_false_num = c(all_phenos_false_num, phenos_false_num)
  }
  return (list(all_phecodes, all_phenos_true_num, all_phenos_false_num))
}

predixcan_list_command = function(){
  phecode_list <- paste0("\"", paste(phecode_with_num, collapse = "\" \""), "\"")
  cat("phecode_list=(", phecode_list, ")", sep = "")
  cat("sample_name=\"NF1_disease\"")
  cat("for phecode in ${phecode_list[@]}; do")
  cat(paste0("python3 ~/airr/othercode/MetaXcan/software/PrediXcanAssociation.py ",
              "--expression_file ~/airr/",stratify_data_fn,"/expr_${sample_name}.txt ",
              "--input_phenos_file ~/airr/",stratify_data_fn,"/X${phecode}_${sample_name}.tsv ",
              "--input_phenos_column X${phecode} ",
              "--covariates_file ~/airr/",stratify_data_fn,"/cov_${sample_name}.txt ",
              "--covariates isMale PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 age ",
              "--output ~/airr/",stratify_result_fn,"/assoc_${phecode}_${sample_name}.tsv ", 
              "--mode logistic"))
  cat("done")
}

pathogenic_variants_score = function(rng, ukb_gene_var_fn, score_coef_fn, vep_fn){
  # Create score coef file for score calculation
  
  # Get pathogenic variants ID from ClinVar
  fl = VcfFile("interpretation/clinvar/clinvar_20230604.vcf.gz")
  clinvar_NF1 = VariantAnnotation::readVcf(fl, "hg38", param=rng)
  
  clinvar_NF1_info = VariantAnnotation::info(clinvar_NF1)
  clinvar_NF1_info = as.data.frame(clinvar_NF1_info)
  print(base::table(unlist(clinvar_NF1_info$CLNSIG)))
  print(base::table(unlist(clinvar_NF1_info$CLNVC)))

  clinvar_NF1_pathogenic = clinvar_NF1[
    (clinvar_NF1_info$CLNSIG %in% c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic")) 
    # & (clinvar_NF1_info$CLNVC == "single_nucleotide_variant")
    ]
  pathogenic_vid = paste(
    as.character(seqnames(clinvar_NF1_pathogenic)), 
    as.integer(start(clinvar_NF1_pathogenic)), 
    as.character(ref(clinvar_NF1_pathogenic)), 
    sapply(alt(clinvar_NF1_pathogenic), as.character), 
    sep="_"
  )

  # Get pathogenic variants ID from VEP
  vep = read_tsv(vep_fn, comment = "##", name_repair = "universal")
  vep = vep %>% select(.Uploaded_variation, Consequence, IMPACT, SYMBOL)
  pathogenic_vep = vep %>% filter(IMPACT == "HIGH")
  print(table(pathogenic_vep$Consequence))
  print(table(pathogenic_vep$SYMBOL))
  pathogenic_vep = unique(pathogenic_vep$.Uploaded_variation)
  
  # Get UKB variants ID
  ukb_NF1_var = read_tsv(ukb_gene_var_fn, comment = "##")
  ukb_NF1_vid = ukb_NF1_var$ID # base::paste(ukb_NF1_var[["#CHROM"]], ukb_NF1_var$POS, ukb_NF1_var$REF, ukb_NF1_var$ALT, sep="_")
  print(paste("Number of variants overlaping with ClinVar P/LP variants:", length(intersect(pathogenic_vid, ukb_NF1_vid))))
  print(paste("Number of variants overlaping with VEP LoF variants:", length(intersect(pathogenic_vep, ukb_NF1_vid))))


  intersect_ukb_name = intersect(union(pathogenic_vid, pathogenic_vep), ukb_NF1_vid)

  score_coef = data.frame(
    variant = intersect_ukb_name,
    allele = sapply(strsplit(intersect_ukb_name, "_"), function(x) x[[4]]),
    score = rep(1, length(intersect_ukb_name))
  )
  print(head(score_coef))
  print(dim(score_coef))
  write_tsv(score_coef, score_coef_fn)

}

generate_vcf_from_snp_list <- function(snp_list, vcf_filename){
  snp_list = data.frame(ID = snp_list, stringsAsFactors = FALSE)
  snp_list <- separate(snp_list, ID, into = c("chromosome", "position", "ref_allele", "alt_allele"), sep = "_", remove = FALSE)  %>%
    mutate(QUAL = ".", FILTER = ".", INFO = ".")
  vcf_sorted <- snp_list %>% arrange(chromosome, position) %>%
    dplyr::select(chromosome, position, ID, ref_allele, alt_allele, QUAL, FILTER, INFO)
  write_tsv(vcf_sorted, vcf_filename, col_names = FALSE)
}
