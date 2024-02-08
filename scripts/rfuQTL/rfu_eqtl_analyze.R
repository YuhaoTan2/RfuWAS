# Reference genome: hg38
setwd("~/airr")

library(regioneR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(karyoploteR)
library(MatrixEQTL)
library(tidyverse)
library(vroom)
library(ggplot2)
library(ggpubr)
library(ggsci)
source("scripts/rfu/utils.R")
rename = dplyr::rename


##### qq plot ------------------------------------
prefix = "results/matrixeqtl/"
output_prefix = "tcrnum_residuals"
output_file_name = paste(prefix, "matrixeqtl.", output_prefix, ".tsv", sep="")
rds_file_name = paste(prefix, "matrixeqtl.", output_prefix, ".rds", sep="") 

me = readRDS(rds_file_name)
svg(filename = "figs/qqplot.svg", width = 4, height = 4)
plot(me, main = "rfuQTL QQ-plot")
dev.off()


# eqtl analysis
eqtls = as_tibble(read.table(output_file_name, sep="\t", header = 1))
dim(eqtls)
tail(sort(table(eqtls$gene)))
length(unique(eqtls$gene))
#hist(table(eqtls$gene), breaks=50)
eqtls$CHR = paste("chr", sapply(strsplit(eqtls$SNP,"_"), `[`, 1), sep = "")
eqtls$BP = as.numeric(sapply(strsplit(eqtls$SNP,"_"), `[`, 2))


# Count significant number
# length(unique((eqtls %>% filter(p.value < 0.05/4953/6390032))$SNP))
eqtl.sig = eqtls %>% filter(p.value < 5e-8/4953)
write_tsv(eqtl.sig %>% arrange(p.value) %>% mutate(gene = gene + 1), "tables/rfuqtl.tsv")
cat("Significant RFU-variant:", nrow(eqtl.sig))
table(eqtl.sig$CHR)
cat("Significant variant:", length(unique(eqtl.sig$SNP)))
cat("Significant RFU:", length(unique(eqtl.sig$gene)))

##### MatrixEQTL: Manhattan plot for all eqtl ------------------------------------

# Function to perform thinning on the data
thin_data <- function(data, thin_pos_places = 0, thin_logp_places = 0) {
  data$rounded_pos <- round(data$BP / 10^6, digits = thin_pos_places) * 10^6
  data$rounded_pval <- round(-log10(data$p.value), digits = thin_logp_places)
  thinned_data <- data[!duplicated(data[, c("CHR", "rounded_pos", "rounded_pval")]), ]
  return(thinned_data)
}

# Apply thinning to your data
eqtls_whole = as_tibble(read.table(paste(prefix, "matrixeqtl.", "tcrnum_residuals_1e5", ".tsv", sep=""), sep="\t", header = 1))
eqtls_whole$CHR = paste("chr", sapply(strsplit(eqtls_whole$SNP,"_"), `[`, 1), sep = "")
eqtls_whole$BP = as.numeric(sapply(strsplit(eqtls_whole$SNP,"_"), `[`, 2))
eqtls_thinned <- thin_data(eqtls_whole) %>% arrange(p.value)
dim(eqtls_thinned)

snps = makeGRangesFromDataFrame(eqtls_thinned, seqnames.field = "CHR", start.field = "BP",
                                end.field = "BP", ignore.strand = T)
snps$pval = eqtls_thinned$p.value
snps = setNames(snps, eqtls_thinned$SNP)

ymax = ceiling(-log10(min(snps$pval))) + 15

svg(filename = "figs/manhattan.svg", width = 8, height = 4)
kp <- plotKaryotype(genome='hg38', plot.type=4, labels.plotter = NULL, ideogram.plotter=NULL)
kpAddChromosomeNames(kp, chr.names = c(1:12, "", 14, "", 16, "", 18, "", 20, "", 22, "X", "Y"), cex=1, yoffset = 23)

kp <- kpPlotManhattan(kp, data=snps, ymax=ymax, 
                      genomewideline = -log10(5e-8/4953), points.col = "2blues",
                      suggestiveline = 0, suggestive.lwd = 0)
kpAxis(kp, ymin = 0, ymax=ymax, tick.pos = c(0, 10, 20, 30, 40, 50, 60, 70))
kpAbline(kp, h = 0)

snps <- kp$latest.plot$computed.values$data
#Get the names of the top SNP per chr
top.snps <- tapply(seq_along(snps), seqnames(snps), function(x) {
  in.chr <- snps[x]
  top.snp <- in.chr[which.max(in.chr$y)]
  return(names(top.snp))
})

# For normalize.cov7, normalize.cov9, residuals
chr.to.plot = c("chr7", "chr6")
#labels.to.plot = c("TRBV12-5", "HLA-DPA1")
labels.to.plot = c("TRB", "HLA")

kpPlotMarkers(kp, data=snps[top.snps[chr.to.plot]], labels=labels.to.plot, 
              srt=45, y=8, 
              ymax=ymax, line.color="red", r0 = 0.85)
kpSegments(kp, data=snps[top.snps[chr.to.plot]], y0=snps[top.snps[chr.to.plot]]$y, 
           y1=0.85*ymax, ymax=ymax, col="red")

dev.off()



## Regulation region ----------------------
### cCREs
# instructions: https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1707106066_B7K39Vcy6DbEnBvJvliThUwQikM9&db=hg38&c=chr7&g=encodeCcreCombined
# TRB with 1Mb region hg38: 7:141299011-143813287
cal_enrich_pvalue = function(occur_count, groupby_name, total_eqtl, total_noneqtl){
  ccreLabel_occur = occur_count %>% group_by(!!sym(groupby_name)) %>% 
    summarise(eqtl_count = sum(eqtl_count), all_count = sum(all_count), noneqtl_count = sum(noneqtl_count))
  
  ccreLabel_enrich <- ccreLabel_occur %>% 
    rowwise() %>% 
    mutate(
      p_value = fisher.test(matrix(
        c(eqtl_count, noneqtl_count, 
          total_eqtl - eqtl_count, total_noneqtl - noneqtl_count), 
        ncol = 2, byrow = TRUE),  alternative = "greater")$p.value
    )
  ccreLabel_enrich %>% arrange(p_value)
}

custom_format <- function(x) {
  # Remove the first 0 in 3e-05
  sci_format = x
  if (grepl("e", sci_format)) {
    # Extract the coefficient and the exponent
    parts <- unlist(strsplit(sci_format, "e-"))
    coeff <- parts[1]
    exponent <- parts[2]
    
    # If the exponent has a leading 0, remove it
    if (substr(exponent, 1, 1) == "0") {
      exponent <- substr(exponent, 2, 2)
    }
    return(paste0(coeff, "e-", exponent))
  } else {
    return(sci_format)
  }
}

compare_pvalue = function(data, xlabel, xname, total_eqtl, total_noneqtl, filter = T){
  if (filter) data = data %>% filter(p_value < 0.01)
  data <- data %>% 
    mutate(eqtl_count = eqtl_count / total_eqtl, noneqtl_count = noneqtl_count / total_noneqtl) %>%
    mutate(significance = case_when(
      p_value < 1e-4 ~ '***',
      p_value < 1e-3  ~ '**',
      p_value < 1e-2  ~ '*',
      #p_value < pvalue.thre ~ custom_format(format(p_value, digits=1, scientific=T)),
      TRUE            ~ ''
    ),
    position = max(eqtl_count, noneqtl_count))
  
  data <- data %>%
    arrange(p_value)
  
  data[[xname]] <- factor(data[[xname]], levels = unique(data[[xname]]))
  
  data_long <- data %>%
    gather(key = "type", value = "count", eqtl_count, noneqtl_count) %>%
    mutate(type = recode(type, 
                         "eqtl_count" = "rfuQTL", 
                         "noneqtl_count" = "non-rfuQTL"))

  ggplot(data_long , aes(x = !!sym(xname), y = count, fill = type)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    geom_text(data = filter(data_long, type == "rfuQTL"), aes(x = !!sym(xname), y = position, label = significance), vjust = -0.5) +
    labs(title = "",
         x = xlabel,
         y = "Overlap fraction",
         fill = "Variants") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.8)), 
          axis.text.y = element_text(size = rel(0.8)), 
          axis.title = element_text(size = rel(0.8)), 
          legend.title = element_text(size = rel(0.8)), 
          legend.text = element_text(size = rel(0.8))) + 
    coord_cartesian(ylim = c(min(data_long$count, na.rm = TRUE), max(data_long$count, na.rm = TRUE)*1.1)) + scale_fill_d3()
}


ccre = read_tsv("interpretation/annotation/trb_encodeccrecombined.bed", col_names = F)
colnames(ccre) = c("chrom", "chromStart", "chromEnd", "name", "score",
                   "strand", "thickStart", "thickEnd", "reserved", "ccre",
                    "encodeLabel", "zScore", "ucscLabel", "accessionLabel", "description")


trb_snps = read_tsv("data/gwas/sel_trb.pvar", comment = "##")
trb_snps = trb_snps %>% select(`#CHROM`, POS, ID, REF, ALT) %>% distinct()
hla_snps = read_tsv("data/gwas/sel_hla.pvar", comment = "##")
hla_snps = hla_snps %>% select(`#CHROM`, POS, ID, REF, ALT) %>% distinct()
trb_hla_snps = rbind(trb_snps, hla_snps)

eqtl_snps = eqtl.sig %>% filter(CHR == "chr7") %>% select(SNP, CHR, BP) %>% distinct()

ccre_occur <- ccre %>%
  rowwise() %>%
  mutate(eqtl_count = sum(eqtl_snps$BP >= (chromStart + 1) & eqtl_snps$BP <= (chromEnd + 1), na.rm = TRUE),
         all_count = sum(trb_snps$POS >= (chromStart + 1) & trb_snps$POS <= (chromEnd + 1), na.rm = TRUE),
         noneqtl_count = all_count - eqtl_count,
         width = chromEnd - chromStart + 1) %>%
  ungroup() %>%
  select(chrom, chromStart, chromEnd, name, ccre, encodeLabel, eqtl_count, all_count, noneqtl_count, width)
ccre_occur


total_eqtl <- nrow(eqtl_snps)
total_noneqtl <- nrow(trb_snps)
encodeLabel_pvalue = cal_enrich_pvalue(ccre_occur, "encodeLabel", total_eqtl, total_noneqtl)

svg("figs/encode_ccre.svg",  width = 3.5, height = 3)
compare_pvalue(encodeLabel_pvalue, "ENCODE cCREs", "encodeLabel", total_eqtl, total_noneqtl, filter=F)
dev.off()

## tfbs
tfbs = read_tsv("interpretation/annotation/encRegTfbsClusteredWithCells.hg38.bed", col_names = F)
colnames(tfbs) = c("chrom", "chromStart", "chromEnd", "name", "score", "celltype")
tfbs_trb = tfbs %>% filter(chrom == "chr7", chromStart + 1 >= 141299011, chromEnd + 1 <= 143813287)


tfbs_occur <- tfbs_trb %>%
  rowwise() %>%
  mutate(eqtl_count = sum(eqtl_snps$BP >= (chromStart + 1) & eqtl_snps$BP <= (chromEnd + 1), na.rm = TRUE),
         all_count = sum(trb_snps$POS >= (chromStart + 1) & trb_snps$POS <= (chromEnd + 1), na.rm = TRUE),
         noneqtl_count = all_count - eqtl_count) %>%
  ungroup()

tfbs_pvalue = cal_enrich_pvalue(tfbs_occur, "name", total_eqtl, total_noneqtl)
tfbs_pvalue %>% filter(p_value < 0.05) %>% print(n = Inf)
tfbs_pvalue %>% filter(p_value < 0.05 / nrow(tfbs_pvalue)) %>% print(n = Inf)

svg("figs/encode_tfbs.svg",  width = 5, height = 3)
compare_pvalue(tfbs_pvalue, "ENCODE Transcription Factors", "name", total_eqtl, total_noneqtl, 0.05 / nrow(tfbs_pvalue))
dev.off()


##### HLA association --------------
# Load cd4_cd8_m
cd4_cd8 = read_tsv("results/rfu/cd4_cd8_ms_zumla.tsv")
cd4_cd8_m = cd4_cd8 %>% rename(gene = RFU) %>% select(enrichment, gene, median_ratio.ms, median_ratio.zumla, p.value.ms, p.value.zumla) %>% mutate(median_ratio = (median_ratio.ms + median_ratio.zumla) / 2) %>% filter(p.value.ms < 0.05, p.value.zumla < 0.05)
convert_to_sci = function(pvalue){
  # Split into base and exponent
  parts <- strsplit(sprintf("%.1e", pvalue), "e")[[1]]
  base <- parts[1]
  exponent <- as.numeric(parts[2])
  
  # Create an expression for the label
  bquote(p == .(base) %*% 10^.(exponent))
}

## For all comparision, we need to filter based on AF of TEDDY and Soumya
teddy.freq = read_tsv("data/TEDDY/HLA_imputation/chr6.R20.7.dose.sel.afreq")
teddy.freq = teddy.freq %>% rename(SNP = ID, freq.teddy = ALT_FREQS)
teddy.freq = teddy.freq %>% select(SNP, freq.teddy)
hla_genotype_t = read_tsv("data/elife_discovery/hla_matrixeqtl.traw")
hla_af = hla_genotype_t %>%
  rowwise() %>%
  mutate(af = mean(c_across(starts_with("P")), na.rm = T)) %>%
  select(SNP, af)
hla_af = hla_af %>% mutate(SNP = str_replace(SNP, "^HLA_([A-Z0-9]+)_([0-9]{2})_([0-9]{2,3})$", "HLA_\\1*\\2:\\3"))
teddy.soumya.af = inner_join(teddy.freq, hla_af)

hla.teddyhf = teddy.soumya.af %>% filter(freq.teddy > 0.05) %>% pull(SNP)
hla.soumyahf = teddy.soumya.af %>% filter(af > 0.05) %>% pull(SNP)
hla.highfreq = teddy.soumya.af %>% filter(freq.teddy > 0.05, af > 0.05) %>% pull(SNP)
hla.consistence = teddy.soumya.af %>% filter(abs(freq.teddy - af) < 0.2) %>% pull(SNP)

# TEDDY HLA association
teddy_hla_eqtls = vroom(paste(prefix, "matrixeqtl.", "hla_tcrnum_residuals", ".tsv", sep=""), .name_repair = "universal")
teddy_hla_eqtls = teddy_hla_eqtls %>% filter(SNP %in% hla.teddyhf)
eqtls_rfu_num = length(unique(teddy_hla_eqtls$gene))
eqtls_allele_num = length(unique(teddy_hla_eqtls$SNP))
p.teddy.thre = 0.05 / eqtls_rfu_num / eqtls_allele_num
teddy_hla_eqtls.sig = teddy_hla_eqtls %>% filter(p.value < 0.05 / eqtls_rfu_num / eqtls_allele_num)
teddy_hla_eqtls.sig %>% print(n = Inf)
write_tsv(teddy_hla_eqtls.sig %>% arrange(p.value) %>% mutate(gene = gene + 1), "tables/teddy_hla_qtls.tsv")
teddy_hla_eqtls.sig %>% rename(RFU = gene) %>% mutate(RFU = RFU + 1) %>% left_join(cd4_cd8, by = "RFU") %>% select(SNP, RFU, p.value, beta, enrichment) %>% arrange(p.value) %>% print(n = Inf)

## Compare with CD4 and CD8 enrichment
teddy_mhc2_eqtls = teddy_hla_eqtls %>% 
  filter(grepl(paste0("^", c("HLA_DQA1", "HLA_DQB1", "HLA_DRB1", "HLA_DPA1", "HLA_DPB1"), collapse = "|"), SNP)) %>% 
  group_by(gene) %>% 
  slice_min(p.value) %>% 
  ungroup() %>% 
  select(gene, beta, t.stat, p.value) %>% 
  mutate(gene = gene + 1)
teddy_mhc1_eqtls = teddy_hla_eqtls %>% 
  filter(grepl(paste0("^", c("HLA_A", "HLA_B", "HLA_C"), collapse = "|"), SNP)) %>%
  group_by(gene) %>% 
  slice_min(p.value) %>% 
  ungroup() %>% 
  select(gene, beta, t.stat, p.value) %>% 
  mutate(gene = gene + 1)
teddy_mhc_eqtls = inner_join(teddy_mhc1_eqtls, teddy_mhc2_eqtls, by = "gene", suffix = c(".mhc1", ".mhc2")) %>% 
  mutate(gene = gene + 1, mhc2_mhc1_ratio = abs(beta.mhc2 / beta.mhc1)) %>% 
  left_join(cd4_cd8_m, by = "gene") %>% filter(median_ratio != 1)# %>%
  # filter(p.value.zumla < 0.05, median_ratio.zumla != 1)

# Scatter plots
cor_results <- cor.test(abs(teddy_mhc_eqtls$mhc2_mhc1_ratio), teddy_mhc_eqtls$median_ratio, method = 's')
cor_val <- cor_results$estimate
p_val <- cor_results$p.value

cor_val
p_val

ggplot(teddy_mhc_eqtls, aes(x = mhc2_mhc1_ratio, y = (median_ratio))) +
  geom_point(size = 0.5) +
  labs(
    x = "MHC-II/MHC-I effect size ratio",
    y = "CD4/CD8 ratio",
    title = ""
  ) +
  geom_smooth(method = 'lm') +
  annotate("text", x = 3.5, y = 5, 
           label = sprintf("r = %.2f", cor_val), colour = "red") +
  annotate("text", x = 3.5, y = 4.5, 
           label = convert_to_sci(p_val), colour = "red") +
  theme_minimal()  # Optional: Use a minimal theme
ggsave("figs/compare_teddy_mhc_cd4_ratio.svg",  width = 3, height = 3)




##### Soumya HLA association -------------------
soumya_hla_eqtls = vroom(paste(prefix, "matrixeqtl.", "hla_soumya_2e4_residuals", ".tsv", sep=""), .name_repair = "universal")
soumya_hla_eqtls = soumya_hla_eqtls %>% 
  mutate(SNP = str_replace(SNP, "^HLA_([A-Z0-9]+)_([0-9]{2})_([0-9]{2,3})$", "HLA_\\1*\\2:\\3")) %>% 
  filter(SNP %in% hla.soumyahf)
eqtls_rfu_num = length(unique(soumya_hla_eqtls$gene))
eqtls_allele_num = length(unique(soumya_hla_eqtls$SNP))
soumya_hla_eqtls.sig = soumya_hla_eqtls %>% filter(p.value < 0.05 / eqtls_rfu_num / eqtls_allele_num)
write_tsv(soumya_hla_eqtls.sig %>% arrange(p.value) %>% mutate(gene = gene + 1), "tables/soumya_hla_qtls.tsv")
soumya_hla_eqtls.sig %>% rename(RFU = gene) %>% mutate(RFU = RFU + 1) %>% left_join(cd4_cd8, by = "RFU") %>% select(SNP, RFU, p.value, beta, enrichment) %>% arrange(p.value) %>% print(n = Inf)

## Compare TEDDY and Soumya HLA association
hla_eqtls = inner_join(teddy_hla_eqtls, soumya_hla_eqtls, by=c("SNP", "gene"), suffix = c(".teddy", ".soumya"))
hla_eqtls_sel <- hla_eqtls %>%
  filter(SNP %in% hla.highfreq ) %>%
  #filter(grepl(paste0("^", c("HLA_DQA1", "HLA_DQB1", "HLA_DRB1", "HLA_DPA1", "HLA_DPB1"), collapse = "|"), SNP)) %>% 
  #filter(grepl(paste0("^", c("HLA_A", "HLA_B", "HLA_C"), collapse = "|"), SNP)) %>% 
  group_by(SNP) %>%
  #filter(p.value.teddy < p.teddy.thre) %>%
  slice_min(p.value.teddy, n = 5) %>%
  ungroup() 

cor_results <- cor.test(hla_eqtls_sel$beta.teddy, hla_eqtls_sel$beta.soumya)
cor_val <- cor_results$estimate
p_val <- cor_results$p.value

cor_val
p_val

ggplot(hla_eqtls_sel, aes(x = beta.teddy, y = beta.soumya)) +
  geom_point(size = 0.5) +
  labs(
    x = "Effect size on training set",
    y = "Effect size on test set",
    title = ""
  ) +
  geom_smooth(method = 'lm') +
  annotate("text", x = 0.3, y = -0.4, 
           label = sprintf("r = %.2f", cor_val), colour = "red") +
  annotate("text", x = 0.3, y = -0.55, 
           label = convert_to_sci(p_val), colour = "red") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figs/compare_teddy_soumya_hla_beta.svg",  width = 3, height = 3)
