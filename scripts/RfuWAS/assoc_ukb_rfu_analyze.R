# Significant results are saved to results/assoc_pheno/ukb_hg19_phenome_sig.tsv

setwd("~/airr")
library(tidyverse)
library(readr)
library(vroom)
library(ggplot2)
library(ggsci)
library(GWASTools)
"%&%" = function(a,b) paste0(a,b)
source("scripts/rfu/utils.R")
source("scripts/ukb/ukb_utils.R")
source("scripts/ukb/pheno_utils.R")

### Select significant results
collect_assoc = function(phenome.des.fn = "data/ukb/qc/UKB_PHENOME_DESCRIPTION_QC.txt", 
                        dataset.name = "ukbflip_hg19"){
  phenome.des = read_tsv(phenome.des.fn)
  
  output_dir = "results/assoc_pheno/"%&% dataset.name %&% "/"

  assoc.fn = paste(output_dir,"assoc_X", phenome.des$phecode[1],".tsv", sep = "")
  assoc_df = read_tsv(assoc.fn)
  hg19.rfu = assoc_df$gene
  rfu.num = nrow(assoc_df)
  cat("Number of genes: ", dim(assoc_df)[1])

  phenome.sig = data.frame()
  notfound_file = c()

  for (pheno.name in phenome.des$phecode){
    assoc.fn = paste(output_dir, "assoc_X", pheno.name,"_permuted.tsv", sep = "")
    if (!file.exists(assoc.fn)) {cat(pheno.name," file not found\n")
      notfound_file = c(notfound_file, pheno.name)
      next}
    assoc_df = vroom(assoc.fn, show_col_types = FALSE)
    assoc_df_sig = assoc_df %>% 
      arrange(pvalue)
    if (nrow(assoc_df_sig) > 0) {
      pheno.name.des = phenome.des %>% filter(phecode == pheno.name)%>% select(phecode, description)
      assoc_df_sig = bind_cols(assoc_df_sig, pheno.name.des)
      phenome.sig = bind_rows(phenome.sig, assoc_df_sig)
    }
  }
  # pvalue threshold: 0.05 / rfu.num
  phenome.sig = phenome.sig %>% rename(RFU = gene)
  write_tsv(phenome.sig, paste0("results/assoc_pheno/", dataset.name, "_phenome_permuted_all.tsv"))

  # pvalue threshold: 0.05 / rfu.num / nrow(phenome.des)
  #phenome.sig.p = phenome.sig %>% filter(pvalue < 0.05 / rfu.num / nrow(phenome.des)) %>%
  #  arrange(pvalue)
  #head(phenome.sig.p, 10)
  #phenome.sig.p = phenome.sig.p %>% select(RFU, pvalue, effect, phecode, description)
  #write_tsv(phenome.sig.p, paste0("results/assoc_pheno/", dataset.name, "_phenome_sig.tsv"))
  cat("Not found file:", notfound_file, "\n")
  return ()
}

## create the volcano plot
volcano_plot = function(phenome.sig){
  phenome.sig <- phenome.sig %>% mutate(pvalue = ifelse(pvalue == 0, .Machine$double.xmin, pvalue)) %>%
    # mutate(significant = ifelse(pvalue < 0.05 / rfu.num / nrow(phenome_des) & abs(effect) > 0.1, "Sig", "Not Sig"))
    mutate(significant = ifelse(pvalue < 0.05 / rfu.num, "Sig", "Not Sig"))
  ggplot(phenome.sig, aes(x = effect, y = -log10(pvalue), color = significant)) +
    geom_point(size = 1) +
    scale_color_manual(values = c("Not Sig" = "grey", "Sig" = "red")) +
    geom_hline(yintercept = -log10(0.05 / rfu.num), linetype = "dashed", color = "grey50") +
    labs(x = "effect size", y = "-log10(p-value)")+#, color = "Significant") +
    theme_bw()
}

collect_assoc("data/ukb/qc/UKB_PHENOME_DESCRIPTION_QC.txt", "ukbflip_hg19")


phenome.p.all = read_tsv(paste0("results/assoc_pheno/", "ukbflip_hg19", "_phenome_all.tsv"))
svg("figs/rfuwas_qqplot_ori.svg", width = 3.5, height = 3.5)
qqPlot(phenome.p.all$pvalue, main = "RfuWAS p-values", thinThreshold = -log10(1e-15))
dev.off()

######################################################################
# Start of association analysis, p value threshold 0.05 / rfu.num --------
phenome.p = read_tsv(paste0("results/assoc_pheno/", "ukbflip_hg19", "_phenome.tsv")) # %>% filter(gene %in% rfu_hq)
rfu.num = get_rfu_num() 

phenome_des = vroom("data/ukb/qc/UKB_PHENOME_DESCRIPTION_QC.txt")
phenome_des_notpass = vroom("data/ukb/qc/UKB_PHENOME_DESCRIPTION_NOTPASS.txt")
phenome_des_all = rbind(phenome_des, phenome_des_notpass)
phenome_group_all = phenome_des_all %>% select(phecode, group)
phenome.p = phenome.p %>% left_join(phenome_group_all, by = "phecode")

# Get autoimmune disease
# autoimmune_name.txt is from https://autoimmune.org/disease-information/
autoimmune_name = read_tsv("interpretation/autoimmune_name.txt", col_names = F, show_col_types = FALSE)
# Remove (xxxx) in the disease name
autoimmune_name = (autoimmune_name %>% mutate(X1 = gsub(" \\(.*\\)", "", X1)))$X1
autoimmune_name = c(autoimmune_name, "regional enteritis",
                    "celiac", "Psoriatic arthropathy", "Graves' disease")

# Assign group
phenome.p = phenome.p %>% mutate(group = ifelse(grepl(paste(autoimmune_name, collapse = "|"), description,ignore.case=TRUE), "autoimmune", group))
phenome.p.sig = phenome.p %>% filter(pvalue < 0.05 / rfu.num / nrow(phenome_des))
length(unique(phenome.p.sig$RFU))
length(unique(phenome.p.sig$phecode))
phenome.p.sig %>% mutate(RFU = RFU + 1) %>% select(-status) %>% arrange(pvalue)
sort(table(phenome.p.sig$description))
tail(sort(table(phenome.p.sig$RFU + 1)))
intersect(phenome.p.sig %>% filter(group == "neoplasms") %>% pull(RFU), c(3826, 632, 3725, 4664, 3458))
cd4_cd8 %>% filter(RFU == 4499) %>% select(RFU, enrichment)
treg_enrichment %>% filter(RFU == 4499)
phenome.p.sig %>% filter(RFU == 4498) %>% print(n = Inf)

write_tsv(phenome.p.sig %>% mutate(RFU = RFU + 1) %>% select(-status) %>% arrange(pvalue), "tables/rfuwas.tsv")
phenome_group = phenome_des %>% mutate(group = ifelse(grepl(paste(autoimmune_name, collapse = "|"), description,ignore.case=TRUE), "autoimmune", group))
phenome_group_num = table(phenome_group$group)
phenome_group_num
phenome.sig.group.num = table(phenome.p.sig$group)
model_group_contingency = full_join(as.data.frame(phenome.sig.group.num), 
                                    as.data.frame(phenome_group_num), by = "Var1", suffix = c(".sig", ".all"))
model_group_contingency = model_group_contingency %>% replace_na(list(Freq.sig = 0, Freq.all = 0)) %>%
  mutate(Freq.notsig = Freq.all * rfu.num - Freq.sig)
rownames(model_group_contingency) = model_group_contingency$Var1
model_group_contingency = t(model_group_contingency %>% select(Freq.sig, Freq.notsig))
model_group_test = chisq.test(model_group_contingency)
model_group_test$p.value
model_group_test$residuals

sort(phenome.sig.group.num)
phenome.sig.forplot = phenome.p %>% 
  mutate(pvalue = ifelse(pvalue == 0, .Machine$double.xmin, pvalue)) %>% 
  mutate(group = ifelse(group %in% c("autoimmune", "endocrine/metabolic", "hematopoietic"), group, "Other")) %>% 
  mutate(group = ifelse(group == "endocrine/metabolic", "Metabolic", group)) %>% 
  mutate(group = ifelse(group == "autoimmune", "Autoimmune", group)) %>% 
  mutate(group = ifelse(group == "hematopoietic", "Hematopoietic", group)) %>% 
  arrange(pvalue)
ggplot(phenome.sig.forplot, aes(x = effect, y = -log10(pvalue), color = group)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05 / rfu.num / nrow(phenome_des)), linetype = "dashed", color = "grey50") +
  labs(x = "Effect size", y = "-log10(p-value)", color = "") +
  theme_bw() +
  theme(legend.position = "bottom", legend.spacing.y = unit(0, "cm"), legend.spacing.x = unit(0, "cm")) +
  guides(color = guide_legend(nrow = 2, byrow = F)) +
  scale_fill_manual(values=c("Autoimmune", "Metabolic", "Hematopoietic", "Other"))
ggsave("figs/phenome_volcano_disease.svg", width = 2.5, height = 3.5)

# Cell type ------------------------------
predictable_rfu_ref = get_predictable_rfu() + 1

cd4_cd8 = read_tsv("results/rfu/cd4_cd8_ms_zumla.tsv") %>% filter(RFU %in% predictable_rfu_ref, !is.na(enrichment))
table(cd4_cd8$enrichment)
phenome.p.sig_ref = phenome.p.sig %>% mutate(RFU = RFU + 1) %>% left_join(cd4_cd8, by = "RFU") 
observed_num = phenome.p.sig_ref %>% group_by(phecode, description, group) %>% 
  summarise(CD4_num = sum(enrichment == "CD4-enriched", na.rm = TRUE), 
            CD8_num = sum(enrichment == "CD8-enriched", na.rm = TRUE), 
            total = n()) %>% ungroup() %>% filter(total > 20) %>% arrange(total) %>%
  filter(!phecode %in% c(244, 250.13, 557, 696))
sel_phecode = unique(observed_num$phecode)

expected_num = table(cd4_cd8$enrichment)
observed_num = observed_num %>% mutate(
  CD4_pvalue = 1 - phyper(CD4_num - 1, expected_num["CD4-enriched"], length(predictable_rfu_ref) - expected_num["CD4-enriched"], total),
  CD8_pvalue = 1 - phyper(CD8_num - 1, expected_num["CD8-enriched"], length(predictable_rfu_ref) - expected_num["CD8-enriched"], total))

write_tsv(observed_num, "tables/rfuwas_cd4_cd8.tsv")

long_data <- observed_num %>% mutate(NA_num = total - CD4_num - CD8_num) %>%
  mutate(CD4_num = CD4_num / expected_num[["CD4-enriched"]], 
         CD8_num = CD8_num / expected_num[["CD8-enriched"]], 
         NA_num = NA_num / (length(predictable_rfu_ref) - sum(expected_num))) %>%
  pivot_longer(cols = c(CD4_num, CD8_num, NA_num), 
               names_to = "cell_type", 
               values_to = "count") %>%
  select(-CD4_pvalue, -CD8_pvalue)

# Add star to plot
annotation_data <- observed_num %>%
  pivot_longer(cols = c("CD4_pvalue", "CD8_pvalue"), names_to = "cell_type", values_to = "pvalue") %>%
  mutate(star = case_when(
    pvalue < 0.005 ~ "**",
    pvalue < 0.05 ~ "*",
    TRUE ~ ""
  )) %>% select(phecode, description, pvalue, star, cell_type)

annotation_data$cell_type <- sub("_pvalue", "_num", annotation_data$cell_type)

annotation_data <- inner_join(annotation_data, long_data)

# Create the histogram using ggplot2
ggplot(long_data, aes(x = count, y = as_factor(description), fill = as_factor(cell_type))) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_y_discrete(labels = paste0(rev(c("Celiac disease", "Psoriasis", "Hypothyroidism", "Iron disorders", "Ankylosing spondylitis", "Type 1 diabetes", "Multiple sclerosis", "Thyrotoxicosis")), " (", observed_num$total, ")")
  ) +
  scale_fill_d3(labels = c(
    "CD4_num" = paste0("CD4 (",expected_num[["CD4-enriched"]],")"), 
    "CD8_num" = paste0("CD8 (",expected_num[["CD8-enriched"]], ")"), 
    "NA_num" = paste0("Unknown\n(", length(predictable_rfu_ref) - sum(expected_num), ")"))) +
  labs(x = "Associated fraction", y = "", fill = "Cell Type (n)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) + 
  geom_text(data = annotation_data, aes(x = count, y = description, label = star),
                                                                       position = position_dodge(width = 0.9), vjust = -0.5)
ggsave("figs/rfuwas_cd4_cd8.svg", width = 4, height = 3.5)

# Plot for positive effect
phenome.p.sig_ref = phenome.p.sig %>% mutate(RFU = RFU + 1) %>% left_join(cd4_cd8, by = "RFU") %>% filter(effect > 0)
observed_num = phenome.p.sig_ref %>% group_by(phecode, description, group) %>% 
  summarise(CD4_num = sum(enrichment == "CD4-enriched", na.rm = TRUE), 
            CD8_num = sum(enrichment == "CD8-enriched", na.rm = TRUE), 
            total = n()) %>% ungroup() %>% filter(phecode %in% sel_phecode) %>% arrange(total)

expected_num = table(cd4_cd8$enrichment)
observed_num = observed_num %>% mutate(
  CD4_pvalue = 1 - phyper(CD4_num - 1, expected_num["CD4-enriched"], length(predictable_rfu_ref) - expected_num["CD4-enriched"], total),
  CD8_pvalue = 1 - phyper(CD8_num - 1, expected_num["CD8-enriched"], length(predictable_rfu_ref) - expected_num["CD8-enriched"], total))

write_tsv(observed_num, "tables/rfuwas_cd4_cd8_positive.tsv")

long_data <- observed_num %>% mutate(NA_num = total - CD4_num - CD8_num) %>%
  mutate(CD4_num = CD4_num / expected_num[["CD4-enriched"]], 
         CD8_num = CD8_num / expected_num[["CD8-enriched"]], 
         NA_num = NA_num / (length(predictable_rfu_ref) - sum(expected_num))) %>%
  pivot_longer(cols = c(CD4_num, CD8_num, NA_num), 
               names_to = "cell_type", 
               values_to = "count") %>%
  select(-CD4_pvalue, -CD8_pvalue)

# Add star to plot
annotation_data <- observed_num %>%
  pivot_longer(cols = c("CD4_pvalue", "CD8_pvalue"), names_to = "cell_type", values_to = "pvalue") %>%
  mutate(star = case_when(
    pvalue < 0.025 ~ "*",
    TRUE ~ ""
  )) %>% select(phecode, description, pvalue, star, cell_type)

annotation_data$cell_type <- sub("_pvalue", "_num", annotation_data$cell_type)

annotation_data <- inner_join(annotation_data, long_data)

# Create the histogram using ggplot2
ggplot(long_data, aes(x = count, y = as_factor(description), fill = as_factor(cell_type))) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_y_discrete(labels = paste0(rev(c("Celiac disease", "Hypothyroidism", "Psoriasis", "Type 1 diabetes", "Iron disorders", "Ankylosing spondylitis", "Thyrotoxicosis", "Multiple sclerosis")), " (", observed_num$total, ")")
  ) +
  scale_fill_d3(labels = c(
    "CD4_num" = paste0("CD4 (",expected_num[["CD4-enriched"]],")"), 
    "CD8_num" = paste0("CD8 (",expected_num[["CD8-enriched"]], ")"), 
    "NA_num" = paste0("Unknown\n(", length(predictable_rfu_ref) - sum(expected_num), ")"))) +
  labs(x = "Associated fraction", y = "", fill = "Cell Type (n)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) + 
  geom_text(data = annotation_data, aes(x = count, y = description, label = star),
            position = position_dodge(width = 1), vjust = 0.3)
ggsave("figs/rfuwas_cd4_cd8_positive.svg", width = 4, height = 3.5)

# Treg ----------------------------------------------------------
treg_enrichment = as_tibble(readRDS("results/rfu/treg_enrichment.rds")) %>% 
  filter(RFU %in% predictable_rfu_ref) %>% 
  mutate(adjusted_Friedman_p = p.adjust(Friedman_p, method = "bonferroni")) %>% 
  filter(adjusted_Friedman_p < 0.05) %>% 
  arrange(Enriched_Group, Friedman_p)
table(treg_enrichment$Enriched_Group)
write_tsv(treg_enrichment, "tables/treg_enrichment_filter.tsv")
# START OF TREG VISUAL ----------------------------------------------------------
library(ComplexHeatmap)
library(circlize)

# Visualize Treg annotation
long_data <- treg_enrichment %>% 
  mutate(RFU = row_number(),
    sum_median = CM_median + TN_median + Treg_median + Tscm_median,
         CM = CM_median / sum_median,
         TN = TN_median / sum_median,
         Treg = Treg_median / sum_median,
         Tscm = Tscm_median / sum_median) %>%
  pivot_longer(cols = c(CM, TN, Treg, Tscm), 
               names_to = "CellType", 
               values_to = "MedianValue")

heatmap_data <- long_data %>%
  pivot_wider(names_from = CellType, values_from = MedianValue) %>%
  column_to_rownames(var = "RFU") %>%
  select(CM, TN, Treg, Tscm)

# Create a color annotation based on Enriched_Group
group_colors <- c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF")
names(group_colors) <- unique(treg_enrichment$Enriched_Group)
row_annotation <- HeatmapAnnotation(df = data.frame(Enriched_Group = treg_enrichment$Enriched_Group),
                                    col = list(Enriched_Group = group_colors),
                                    which = "row",
                                    name = "Cell Type",
                                    annotation_label = "Cell Type",
                                    annotation_name_rot = 45)

heatmap_matrix <- as.matrix(heatmap_data)

svg("figs/treg.svg", width = 3.5, height = 3.5)
Heatmap(heatmap_matrix, 
        name = "Fraction",
        cluster_rows = F, 
        cluster_columns = F,
        show_row_dend = F, 
        show_column_dend = TRUE,
        right_annotation = row_annotation,
        show_row_names = F,
        column_names_rot = 45)
dev.off()

# START OF TREG PROPORTION ----------------------------------------------------------
phenome.p.sig_ref = phenome.p.sig %>% mutate(RFU = RFU + 1) %>% left_join(treg_enrichment, by = "RFU") %>% arrange(pvalue) # %>% filter(effect > 0)

# Normalized RFU counts for each cell types
observed_num = phenome.p.sig_ref %>% group_by(phecode, description, group) %>%
  summarise(CM_num = sum(Enriched_Group  == "CM", na.rm = TRUE), 
            TN_num = sum(Enriched_Group  == "TN", na.rm = TRUE), 
            Tscm_num = sum(Enriched_Group  == "Tscm", na.rm = TRUE), 
            Treg_num = sum(Enriched_Group  == "Treg", na.rm = TRUE), 
            total = n()) %>% ungroup() %>% filter(total > 50) %>% arrange(total) %>%
  filter(!phecode %in% c(244, 250.13, 557, 696))

expected_num = table(treg_enrichment$Enriched_Group)
observed_num = observed_num %>% mutate(
  CM_pvalue = 1 - phyper(CM_num - 1, expected_num["CM"], length(predictable_rfu_ref) - expected_num["CM"], total),
  Treg_pvalue = 1 - phyper(Treg_num - 1, expected_num["Treg"], length(predictable_rfu_ref) - expected_num["Treg"], total),
  TN_pvalue = 1 - phyper(TN_num - 1, expected_num["TN"], length(predictable_rfu_ref) - expected_num["TN"], total),
  Tscm_pvalue = 1 - phyper(Tscm_num - 1, expected_num["Tscm"], length(predictable_rfu_ref) - expected_num["Tscm"], total))

write_tsv(observed_num, "tables/rfuwas_treg.tsv")

long_data <- observed_num %>% mutate(NA_num = total - CM_num - TN_num - Treg_num - Tscm_num) %>%
  mutate(Tscm_num = Tscm_num / expected_num[["Tscm"]], CM_num = CM_num / expected_num[["CM"]], TN_num = TN_num / expected_num[["TN"]], Treg_num = Treg_num / expected_num[["Treg"]], NA_num = NA_num / (length(predictable_rfu_ref) - sum(expected_num))) %>%
  pivot_longer(cols = c(Tscm_num, CM_num, TN_num, Treg_num, NA_num), 
               names_to = "cell_type", 
               values_to = "count") %>%
  select(-CM_pvalue, -TN_pvalue, -Treg_pvalue, -Tscm_pvalue)

# Add star to plot
annotation_data <- observed_num %>%
  pivot_longer(cols = c("CM_pvalue", "TN_pvalue", "Treg_pvalue", "Tscm_pvalue"), names_to = "cell_type", values_to = "pvalue") %>%
  mutate(star = case_when(
    pvalue < 0.0125 ~ "*",
    TRUE ~ ""
  )) %>% select(phecode, description, pvalue, star, cell_type)

annotation_data$cell_type <- sub("_pvalue", "_num", annotation_data$cell_type)

annotation_data <- inner_join(annotation_data, long_data)

# Create the histogram using ggplot2
ggplot(long_data, aes(x = count, y = as_factor(description), fill = as_factor(cell_type))) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_y_discrete(labels = paste0(rev(c("Celiac disease", "Psoriasis", "Hypothyroidism", "Iron disorders", "Ankylosing spondylitis", "Type 1 diabetes", "Multiple sclerosis", "Thyrotoxicosis")), " (", observed_num$total, ")")
                   ) +
  scale_fill_d3(labels = c("CM_num" = paste0("CM (",expected_num[["CM"]],")"), "TN_num" = paste0("TN (",expected_num[["TN"]], ")"), "Tscm_num" = paste0("Tscm (", expected_num["Tscm"], ")"), "Treg_num" = paste0("Treg (", expected_num["Treg"], ")"), "NA_num" = paste0("Unknown\n(", length(predictable_rfu_ref) - sum(expected_num), ")"))) +
  labs(x = "Associated fraction", y = "", fill = "Cell Type (n)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) + 
  geom_text(data = annotation_data, aes(x = count, y = description, label = star),
            position = position_dodge(width = 0.5), 
            vjust = 0.3)
ggsave("figs/rfuwas_celltype.svg", width = 4.5, height = 3.5)



# Antigen-specific TCR in McPAS-TCR -----------------------
annotate_mcpas = function(){
  mcpas = read_csv(paste0(fn.prefix, "McPAS-TCR.csv"))
  mcpas.human = mcpas %>% filter(Species == "Human") %>%
    select(CDR3.alpha.aa, CDR3.beta.aa, Category, Pathology, Antigen.protein, MHC, 
           T.Cell.Type, TRAV, TRAJ, TRBV, TRBD, TRBJ, Antigen.identification.method)
  
  mcpas.human.cdr3 = unique(na.omit(mcpas.human$CDR3.beta.aa))
  mcpas.human.rfu = AssignRFUs.rank.dd(mcpas.human.cdr3)
  mcpas.human.rfu.df = data.frame(CDR3.beta.aa = names(mcpas.human.rfu), rfu = unlist(mcpas.human.rfu))
  mcpas.human = mcpas.human %>% left_join(mcpas.human.rfu.df)
  write_tsv(mcpas.human, "data/antigen_tcr/McPAS-TCR_human_rfu.tsv")
}

mcpas.human = read_tsv(paste0("data/antigen_tcr/", "McPAS-TCR_human_rfu.tsv")) %>% filter(Antigen.identification.method %in% c(1, 2.1, 2.2, 2.3, 2.4, 2.5))
predictable_rfu_ref = get_predictable_rfu() + 1

# rfuWAS can be used to prioritize antigen-specific RFU ------------------------
# Example 1: Celiac disease
mcpas.celiac = mcpas.human %>% filter(Pathology == "Celiac disease")
table(mcpas.celiac$Antigen.protein)
celiac_specific_rfu = intersect(unique(mcpas.celiac$rfu), predictable_rfu_ref)

celiac_rfuwas_ref = read_tsv(paste0("results/assoc_pheno/",  "ukbflip_hg19", "/assoc_X", "557.1",".tsv"), show_col_types = F) %>% 
  rename(RFU = gene) %>% mutate(RFU = RFU + 1) %>%
  filter(RFU %in% predictable_rfu_ref) %>%
  mutate(true_label = ifelse(RFU %in% celiac_specific_rfu, 1, 0))  %>% arrange(pvalue) %>% filter(effect > 0)
celiac_rfuwas_allele1_ref = vroom(paste0("",stratify_result_fn,"/assoc_", "557.1", "_", "dq2.5", ".tsv"), show_col_types = F) %>% 
  rename(RFU = gene) %>% mutate(RFU = RFU + 1) %>%
  filter(RFU %in% predictable_rfu_ref) %>%
  mutate(true_label = ifelse(RFU %in% celiac_specific_rfu, 1, 0)) %>% filter(effect > 0) %>% arrange(pvalue)

celiac.validate.all = read_tsv("results/rfuWmetadata/celiac_RFU_lm.tsv")

celiac_tcr = celiac.validate.all %>% filter(effect > 0) %>% group_by(RFU) %>% slice_min(pvalue) %>% ungroup() %>%  filter(RFU %in% predictable_rfu_ref) %>%
  mutate(true_label = ifelse(RFU %in% celiac_specific_rfu, 1, 0)) %>% arrange(pvalue)


# Compare z-score for antigen-specific and non antigen-specific ------------------------
celiac_rfuwas_ref_sig = celiac_rfuwas_ref %>% filter(pvalue < 0.05 / rfu.num / nrow(phenome_des))
wilcox_result <- wilcox.test(
  celiac_rfuwas_ref_sig %>% filter(true_label == 0) %>% pull(zscore),
  celiac_rfuwas_ref_sig %>% filter(true_label == 1) %>% pull(zscore), alternative = "less"
)

p_value <- wilcox_result$p.value

convert_to_sci = function(p_value){
  # Split into base and exponent
  parts <- strsplit(p_value, "e")[[1]]
  base <- parts[1]
  exponent <- as.numeric(parts[2])
  
  # Create an expression for the label
  label_expr <- bquote(p == .(base) %*% 10^.(exponent))
}

ggplot(celiac_rfuwas_ref_sig, aes(x = as_factor(true_label), y = zscore)) +
  geom_boxplot() +
  scale_x_discrete(labels=c('0' = 'Not specific', '1' = 'CD-specific')) +
  labs(
    x = "gdRFUs", 
    y = "RfuWAS z-score", 
    title = ""
  ) +
  theme_minimal() +
  annotate(
    "text", x = 1.5, y = 1.05 * max(celiac_rfuwas_ref$zscore), 
    label = convert_to_sci(format(p_value, digits = 2, scientific = T)), 
    hjust = 0.5, vjust = 0.2
  ) +theme(plot.title = element_text(hjust = 0.5))+
  coord_cartesian(ylim = c(NA, 1.1 * max(celiac_rfuwas_ref$zscore))) 
ggsave("figs/antigen_celiac.svg", width = 3, height = 2.5)

# sum(celiac_rfuwas_allele2_ref$pvalue < 0.05 / rfu.num)
sum(celiac_rfuwas_ref$pvalue < 0.05 / rfu.num / nrow(phenome_des))


# Plot for all ----------------------------------------------
celiac_specific_freq_pvalue = data.frame()
for (top.num in round(c(10, 20, 50, 100, sum(celiac_rfuwas_ref$pvalue < 0.05 / rfu.num / nrow(phenome_des))))){
  rfuwas_freq = sum(celiac_rfuwas_ref[1:top.num,"true_label"] == 1) / top.num
  rfuwas_pvalue = 1 - phyper(sum(celiac_rfuwas_ref[1:top.num,"true_label"] == 1) - 1, 
                             length(celiac_specific_rfu), 
                             length(predictable_rfu_ref) - length(celiac_specific_rfu), 
                             top.num)
  
  
  tcr_freq = sum(celiac_tcr[1:top.num,"true_label"] == 1) / top.num
  tcr_pvalue = 1 - phyper(sum(celiac_tcr[1:top.num,"true_label"] == 1) - 1, 
                                     length(celiac_specific_rfu), 
                                     length(predictable_rfu_ref) - length(celiac_specific_rfu), 
                                     top.num)
  
  celiac_specific_freq_pvalue = bind_rows(
    celiac_specific_freq_pvalue, bind_rows(
      data.frame(top_num = top.num, method = "RfuWAS", freq = rfuwas_freq, pvalue = rfuwas_pvalue),
      data.frame(top_num = top.num, method = "TCRseq", freq = tcr_freq, pvalue = tcr_pvalue),
    )
  )
}
celiac_specific_freq_pvalue
random_freq = length(celiac_specific_rfu) / length(predictable_rfu_ref)

ggplot(celiac_specific_freq_pvalue, aes(x = factor(top_num), y = freq, fill = method)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(x = "Top N significant gdRFUs", y = "Proportion CD-specific", fill = "", title = "Celiac disease") +
  geom_hline(yintercept = random_freq, linetype = "dashed", color = "#E53528") +
  geom_text(
    aes(label = ifelse(pvalue < 0.05, "*", ""), group = method),
    data = celiac_specific_freq_pvalue,
    position = position_dodge(width = 0.9), 
    vjust = -0.5, 
    size = 3
  ) +
  scale_fill_manual(
    values = c("RfuWAS" = "#55B7E6", "TCRseq" = "#193E8F")
  ) +
  theme(plot.title = element_text(hjust = 0.5)) + theme(
    legend.position = "bottom", 
    legend.direction = "horizontal",
    legend.box = "horizontal" 
  ) +
  coord_cartesian(ylim = c(NA, 1.1 * max(celiac_specific_freq_pvalue$freq))) # Extend the y-axis limits
ggsave("figs/prioritize_celiac.svg", width = 3.5, height = 3)


# Plot for DQ2 ----------------------------------------------
table(mcpas.celiac$MHC)
sum(celiac_rfuwas_allele1_ref$pvalue < 0.05 / rfu.num)
celiac_specific_freq_pvalue_allele1 = data.frame()
for (top.num in round(c(10, 20, 50, 100, sum(celiac_rfuwas_allele1_ref$pvalue < 0.05 / rfu.num)))){
  rfuwas_freq = sum(celiac_rfuwas_ref[1:top.num,"true_label"] == 1) / top.num
  rfuwas_pvalue = 1 - phyper(sum(celiac_rfuwas_ref[1:top.num,"true_label"] == 1) - 1,
                             length(celiac_specific_rfu),
                             length(predictable_rfu_ref) - length(celiac_specific_rfu),
                             top.num)
  
  rfuwas_allele1_freq = sum(celiac_rfuwas_allele1_ref[1:top.num,"true_label"] == 1) / top.num
  rfuwas_allele1_pvalue = 1 - phyper(sum(celiac_rfuwas_allele1_ref[1:top.num,"true_label"] == 1) - 1,
                             length(celiac_specific_rfu),
                             length(predictable_rfu_ref) - length(celiac_specific_rfu),
                             top.num)
  
  celiac_specific_freq_pvalue_allele1 = bind_rows(
    celiac_specific_freq_pvalue_allele1, bind_rows(
      data.frame(top_num = top.num, method = "RfuWAS", freq = rfuwas_freq, pvalue = rfuwas_pvalue),
      data.frame(top_num = top.num, method = "RfuWAS\n(DQ2)", freq = rfuwas_allele1_freq, pvalue = rfuwas_allele1_pvalue),
    )
  )
}
celiac_specific_freq_pvalue_allele1
random_freq = length(celiac_specific_rfu) / length(predictable_rfu_ref)

ggplot(celiac_specific_freq_pvalue_allele1, aes(x = factor(top_num), y = freq, fill = method)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(x = "Top N significant gdRFUs", y = "Proportion CD-specific", fill = "", title = "Celiac disease") +
  geom_hline(yintercept = random_freq, linetype = "dashed", color = "#E53528") +
  geom_text(
    aes(label = ifelse(pvalue < 0.05, "*", ""), group = method),
    data = celiac_specific_freq_pvalue_allele1,
    position = position_dodge(width = 0.9), 
    vjust = -0.5, 
    size = 3
  ) +
  scale_fill_manual(
    values = c("RfuWAS" = "#55B7E6", "RfuWAS\n(DQ2)" = "#F09739")
  ) +
  theme(plot.title = element_text(hjust = 0.5)) + theme(
    legend.position = "bottom", 
    # legend.text = element_text(size = 8), 
    legend.direction = "horizontal",
    legend.box = "horizontal" 
  ) +
  coord_cartesian(ylim = c(NA, 1.1 * max(celiac_specific_freq_pvalue_allele1$freq))) # Extend the y-axis limits
ggsave("figs/prioritize_celiac_allele1.svg", width = 3.5, height = 3)

# Examine CD4 CD8 proportion
celiac_rfuwas_ref %>% filter(pvalue < 0.05 / rfu.num / nrow(phenome_des)) %>% left_join(cd4_cd8 %>% select(RFU, enrichment))



# Example 2: Type 1 diabetes
mcpas.t1d = mcpas.human %>% filter(Pathology == "Diabetes Type 1")
table(mcpas.t1d$Antigen.protein)
write_tsv(rbind(mcpas.celiac, mcpas.t1d), "tables/mcpas_celiac_t1d.tsv")

t1d_specific_rfu = intersect(unique(mcpas.t1d$rfu), predictable_rfu_ref)

t1d_rfuwas_ref = read_tsv(paste0("results/assoc_pheno/",  "ukbflip_hg19", "/assoc_X", "250.1",".tsv"), show_col_types = F) %>% 
  rename(RFU = gene) %>% 
  mutate(RFU = RFU + 1) %>%
  filter(RFU %in% predictable_rfu_ref) %>%
  mutate(true_label = ifelse(RFU %in% t1d_specific_rfu, 1, 0)) %>% 
  filter(effect > 0)
t1d_rfuwas_allele2_ref = vroom(paste0("",stratify_result_fn,"/assoc_", "250.1", "_", "dr4dq8", ".tsv"), show_col_types = F) %>% 
  rename(RFU = gene) %>% mutate(RFU = RFU + 1) %>%
  filter(RFU %in% predictable_rfu_ref) %>%
  mutate(true_label = ifelse(RFU %in% t1d_specific_rfu, 1, 0)) %>% filter(effect > 0)


t1d.validate.all = read_tsv("results/rfuWmetadata/t1d_RFU_lm.tsv")

t1d_tcr = t1d.validate.all %>% # filter(effect > 0)%>% 
  group_by(RFU) %>% slice_min(pvalue) %>% ungroup() %>%  filter(RFU %in% predictable_rfu_ref) %>%
  mutate(true_label = ifelse(RFU %in% t1d_specific_rfu, 1, 0)) %>% arrange(pvalue)

# Compare z-score for antigen-specific and non antigen-specific ------------------------
dim(t1d_rfuwas_ref)
t1d_rfuwas_ref_sig = t1d_rfuwas_ref %>% filter(pvalue < 0.05 / rfu.num / nrow(phenome_des))
dim(t1d_rfuwas_ref_sig)
table(t1d_rfuwas_ref_sig$true_label)
wilcox_result <- wilcox.test(
  t1d_rfuwas_ref_sig %>% filter(true_label == 0) %>% pull(zscore),
  t1d_rfuwas_ref_sig %>% filter(true_label == 1) %>% pull(zscore), alternative = "less"
)

p_value <- wilcox_result$p.value

ggplot(t1d_rfuwas_ref_sig, aes(x = as_factor(true_label), y = zscore)) +
  geom_boxplot() +
  scale_x_discrete(labels=c('0' = 'Not specific', '1' = 'T1D-specific')) +
  labs(
    x = "gdRFUs", 
    y = "RfuWAS z-score", 
    title = ""
  ) +
  theme_minimal() +
  annotate(
    "text", x = 1.5, y = 1.05 * max(t1d_rfuwas_ref$zscore), 
    label = paste("p = ", format(p_value, digits = 2)), 
    hjust = 0.5, vjust = 0
  ) +theme(plot.title = element_text(hjust = 0.5))+
  coord_cartesian(ylim = c(NA, 1.1 * max(t1d_rfuwas_ref$zscore))) 
ggsave("figs/antigen_t1d.svg", width = 3, height = 2.5)

t1d_specific_freq_pvalue = data.frame()
sum(t1d_rfuwas_ref$pvalue < 0.05 / rfu.num / nrow(phenome_des))

sum(t1d_rfuwas_allele2_ref$pvalue < 0.05 / rfu.num / nrow(phenome_des))

# Plot for whole
t1d_specific_freq_pvalue = data.frame()
for (top.num in round(c(5, 10,20, 40, sum(t1d_rfuwas_ref$pvalue < 0.05 / rfu.num / nrow(phenome_des))))){
  rfuwas_freq = sum(t1d_rfuwas_ref[1:top.num,"true_label"] == 1) / top.num
  rfuwas_pvalue = 1 - phyper(sum(t1d_rfuwas_ref[1:top.num,"true_label"] == 1) - 1, 
                             length(t1d_specific_rfu), 
                             length(predictable_rfu_ref) - length(t1d_specific_rfu), 
                             top.num)
  
  
  tcr_freq = sum(t1d_tcr[1:top.num,"true_label"] == 1) / top.num
  tcr_pvalue = 1 - phyper(sum(t1d_tcr[1:top.num,"true_label"] == 1) - 1, 
                          length(t1d_specific_rfu), 
                          length(predictable_rfu_ref) - length(t1d_specific_rfu), 
                          top.num)
  
  t1d_specific_freq_pvalue = bind_rows(
    t1d_specific_freq_pvalue, bind_rows(
      data.frame(top_num = top.num, method = "RfuWAS", freq = rfuwas_freq, pvalue = rfuwas_pvalue),
      data.frame(top_num = top.num, method = "TCRseq", freq = tcr_freq, pvalue = tcr_pvalue),
    )
  )
}
t1d_specific_freq_pvalue
random_freq = length(t1d_specific_rfu) / length(predictable_rfu_ref)

ggplot(t1d_specific_freq_pvalue, aes(x = factor(top_num), y = freq, fill = method)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(x = "Top N significant gdRFUs", y = "Proportion T1D-specific", fill = "", title = "Type 1 diabetes") +
  geom_hline(yintercept = random_freq, linetype = "dashed", color = "#E53528") +
  geom_text(
    aes(label = ifelse(pvalue < 0.05, "*", ""), group = method),
    data = t1d_specific_freq_pvalue,
    position = position_dodge(width = 0.9), 
    vjust = -0.5, 
    size = 3
  ) +
  scale_fill_manual(
    values = c("RfuWAS" = "#55B7E6", "TCRseq" = "#193E8F")
  ) +
  theme(plot.title = element_text(hjust = 0.5)) + theme(
    legend.position = "bottom", 
    # legend.text = element_text(size = 8), 
    legend.direction = "horizontal",
    legend.box = "horizontal" 
  ) +
  coord_cartesian(ylim = c(NA, 1.1 * max(t1d_specific_freq_pvalue$freq))) # Extend the y-axis limits
ggsave("figs/prioritize_t1d.svg", width = 3.5, height = 3)

table(mcpas.t1d$MHC)
sum(t1d_rfuwas_allele2_ref$pvalue < 0.05 / rfu.num)

# Plot for DR4-DQ8
t1d_specific_freq_pvalue_allele2 = data.frame()
for (top.num in round(c(5, 10,20, 40, sum(t1d_rfuwas_allele2_ref$pvalue < 0.05 / rfu.num)))){
  rfuwas_freq = sum(t1d_rfuwas_ref[1:top.num,"true_label"] == 1) / top.num
  rfuwas_pvalue = 1 - phyper(sum(t1d_rfuwas_ref[1:top.num,"true_label"] == 1) - 1, 
                             length(t1d_specific_rfu), 
                             length(predictable_rfu_ref) - length(t1d_specific_rfu), 
                             top.num)
  
  rfuwas_allele2_freq = sum(t1d_rfuwas_allele2_ref[1:top.num,"true_label"] == 1) / top.num
  rfuwas_allele2_pvalue = 1 - phyper(sum(t1d_rfuwas_allele2_ref[1:top.num,"true_label"] == 1) - 1,
                                     length(t1d_specific_rfu),
                                     length(predictable_rfu_ref) - length(t1d_specific_rfu),
                                     top.num)
  
  
  t1d_specific_freq_pvalue_allele2 = bind_rows(
    t1d_specific_freq_pvalue_allele2, bind_rows(
      data.frame(top_num = top.num, method = "RfuWAS", freq = rfuwas_freq, pvalue = rfuwas_pvalue),
      data.frame(top_num = top.num, method = "RfuWAS\n(DR4)", freq = rfuwas_allele2_freq, pvalue = rfuwas_allele2_pvalue),
    )
  )
}
t1d_specific_freq_pvalue
random_freq = length(t1d_specific_rfu) / length(predictable_rfu_ref)

ggplot(t1d_specific_freq_pvalue_allele2, aes(x = factor(top_num), y = freq, fill = method)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(x = "Top N significant gdRFUs", y = "Proportion T1D-specific", fill = "", title = "Type 1 diabetes") +
  geom_hline(yintercept = random_freq, linetype = "dashed", color = "#E53528") +
  geom_text(
    aes(label = ifelse(pvalue < 0.05, "*", ""), group = method),
    data = t1d_specific_freq_pvalue_allele2,
    position = position_dodge(width = 0.9), 
    vjust = -0.5, 
    size = 3
  ) +
  scale_fill_manual(
    values = c("RfuWAS" = "#55B7E6", "RfuWAS\n(DR4)" = "#F09739")
  ) +
  theme(plot.title = element_text(hjust = 0.5)) + theme(
    legend.position = "bottom", 
    # legend.text = element_text(size = 8), 
    legend.direction = "horizontal",
    legend.box = "horizontal" 
  ) +
  coord_cartesian(ylim = c(NA, 1.1 * max(t1d_specific_freq_pvalue_allele2$freq))) # Extend the y-axis limits
ggsave("figs/prioritize_t1d_allele2.svg", width = 3.5, height = 3)


