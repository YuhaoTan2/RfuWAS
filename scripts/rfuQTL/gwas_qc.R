setwd("~/airr")
library(tidyverse)
source("scripts/rfu/utils.R")
TEDDY.pheno = read.table("data/rawdata/TEDDY.pheno.sel.tsv", sep="\t", header=T)
expression_file_name = "data/qctmp/subject_rfu_matrixeqtl.tsv";
subject.RFU = read.table(expression_file_name, sep="\t", header=T, row.names = 1, comment.char = "")

suffix = "tcrnum"

##### 7. PCA --------------------------
pca <- read.table("./data/qctmp/sel_pca.eigenvec", head=TRUE, comment.char = "")[,-1]
eigenval <- scan("./data/qctmp/sel_pca.eigenval")


pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval)*100)
plot(pve, ylab="Percent variance explained", xlab = "Principal component")

# subject source
sample.subject.source = select(TEDDY.pheno, sample_id, subject_source)
sample.subject.source = distinct(sample.subject.source)
rownames(sample.subject.source) = sample.subject.source$sample_id
subject.source = sample.subject.source[pca$IID, "subject_source"]


# predicted population
predicted.pop = read.table("data/qctmp/sel_rel_smps.txt", sep="\t", head=TRUE, comment.char = "")
predicted.pop[predicted.pop$Computed.population=="", "Computed.population"] = "not assigned"

##### 8. TCR number ----------------------------
subject_sample = TEDDY.pheno[,c("sample_id", "subject_id")]
subject_sample = subject_sample[!duplicated(subject_sample),]
rownames(subject_sample) = subject_sample$sample_id
# Select 725 subjects
subject_ids = subject_sample[colnames(subject.RFU), "subject_id"]
length(subject_ids)
tcr.num = cal_tcr_num("data/TEDDY/TEDDY_airr_subject/", paste0("TEDDY_sub_", subject_ids))
#tcr.num = cal_tcr_num("data/TEDDY/TEDDY_airr_subject_trunc1e4/", paste0("TEDDY_sub_", subject_ids))
max(tcr.num)
min(tcr.num)
sum(tcr.num)
mean(tcr.num)
sd(tcr.num)
too.many.tcr = c()
too.few.tcr = colnames(subject.RFU)[tcr.num < 2500]

## Filter and prepare for matrixEQTL
subject.RFU.f = subject.RFU[, !colnames(subject.RFU) %in% c(too.many.tcr, too.few.tcr)]
dim(subject.RFU.f)
write.table(subject.RFU.f, paste0("data/gwas/subject_rfu_matrixeqtl_",suffix,".tsv"), sep='\t', quote=F)

sample.f = data.frame("#FID" = colnames(subject.RFU.f), "IID" = colnames(subject.RFU.f), check.names = F)
write.table(sample.f, paste0("data/gwas/",suffix,"_sample.txt"), sep='\t', quote=F, row.names = F)

# PCA of subjects in subject.RFU
subject.pca <- prcomp(t(subject.RFU.f), scale = TRUE)
ggplot(as.data.frame(subject.pca$x), aes(PC1, PC2)) + geom_point(size = 1) + 
  coord_equal() + theme_light() + ggtitle("RFU PCs after QC") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/rfu_pc_afterqc.svg", width = 3.5, height = 3.5)

subject.pca <- prcomp(t(subject.RFU), scale = TRUE)
ggplot(as.data.frame(subject.pca$x), aes(PC1, PC2)) + geom_point(size = 1) +
 coord_equal() + theme_light() + ggtitle("RFU PCs before QC") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/rfu_pc_beforeqc.svg", width = 3.5, height = 3.5)


##### 9. PEER ----------------------------
## Known covariates
# Genotype PCs
rownames(pca) = pca$IID
known.var = pca[colnames(subject.RFU.f),paste("PC", 1:10, sep="")]

# sex
subject_sex = TEDDY.pheno[,c("sample_id", "sex")]
subject_sex = subject_sex[!duplicated(subject_sex),]
rownames(subject_sex) = subject_sex$sample_id
subject_sex$sex = as.numeric(as.factor(subject_sex$sex))
sex.f = subject_sex[colnames(subject.RFU.f),"sex"] - 1
known.var$sex = sex.f

# ia_case_endptage
subject_ia_age = TEDDY.pheno[,c("sample_id", "ia_case_endptage")]
subject_ia_age = subject_ia_age[!duplicated(subject_ia_age),]
rownames(subject_ia_age) = subject_ia_age$sample_id
ia.age.f = subject_ia_age[colnames(subject.RFU.f),"ia_case_endptage"]
known.var$ia.endptage = ia.age.f

# ia_outcome
subject_ia_outcome = TEDDY.pheno[,c("sample_id", "ia_outcome")]
subject_ia_outcome = subject_ia_outcome[!duplicated(subject_ia_outcome),]
rownames(subject_ia_outcome) = subject_ia_outcome$sample_id
subject_ia_outcome$sex = as.numeric(as.factor(subject_ia_outcome$ia_outcome))
ia.outcome.f = subject_ia_outcome[colnames(subject.RFU.f),"ia_outcome"]
known.var$ia.outcome = ia.outcome.f

# TCR number
subject_ids.f = subject_sample[colnames(subject.RFU.f), "subject_id"]
tcr.num.f = cal_tcr_num("data/TEDDY/TEDDY_airr_subject/", paste0("TEDDY_sub_", subject_ids.f))
sum(tcr.num.f)
mean(tcr.num.f)
sd(tcr.num.f)
known.var$tcr.num = tcr.num.f

# subject source
subject.source.f = as.numeric(as.factor(sample.subject.source[colnames(subject.RFU.f), "subject_source"]))
known.var$source = subject.source.f
known.var$IID = rownames(known.var)
write_tsv(known.var, "data/gwas/known_var.tsv")
cov.for.peer = known.var[, c("PC1", "PC2","PC3", "PC4", "PC5", "sex")]
write.table(cov.for.peer, paste0("data/gwas/subject_rfu_peer_",suffix,".tab"), sep='\t', quote=F, row.names=F, col.names = F)
dim(cov.for.peer)

# correlation between known covariates
GGally::ggcorr(known.var, label = TRUE, label_size = 3, size = 3, layout.exp = 1)


############################
### Normalize and PEER before running code below!
############################
# Calculate pearson and heatmap
# residuals.csv, NxG matrix
# the inferred factors (X.csv, NxK)
# the weights of each factor for every gene (W.csv, GxK)
# the inverse variance of the weights (Alpha.csv, Kx1)
peer.dir = paste0("data/gwas/rfu_normalize_peer_", suffix)
peer_X = read.csv(paste(peer.dir, "/X.csv", sep=""), header=F)

# pearson with imputed covariates
peer_pca_cor = cor(t(peer_X[-1:-7,]), known.var, method = 'p', use="complete.obs")
rownames(peer_pca_cor) = as.numeric(rownames(peer_pca_cor)) - 7
peer_pca_cor2 = as.data.frame(peer_pca_cor) %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
ggplot(as.data.frame(peer_pca_cor2), aes(x = rowname, y = colname, fill = value)) + 
  geom_tile() + xlab("PEER factors") +
  ylab("known covariates") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                       midpoint = 0, 
                       guide = guide_colorbar(title = "PearsonR")) +
  theme_minimal() +
  coord_equal()
# Examine the relationship between PEER10 and tcr.num
plot(known.var$tcr.num, peer_X[8,], xlab="tcr.num", ylab="PEER2")


library(tidyverse)
peer_X_t = as.data.frame(t(peer_X))
ggplot(peer_X_t, aes(V8, V9)) + geom_point(size = 2) + 
  coord_equal() + theme_light()
ggplot(peer_X_t, aes(V11, V12)) + geom_point(size = 2) + 
  coord_equal() + theme_light()
ggplot(peer_X_t, aes(V13, V14)) + geom_point(size = 2) + 
  coord_equal() + theme_light()
ggplot(peer_X_t, aes(V15, V16)) + geom_point(size = 2) + 
  coord_equal() + theme_light()
ggplot(peer_X_t, aes(V17, V18)) + geom_point(size = 2) + 
  coord_equal() + theme_light()
Alpha = read.csv(paste(peer.dir, "/Alpha.csv", sep=""), header=F)$V1
plot(0:10, 1.0 / Alpha[-1:-6],xlab="Intercept and PEER Factors", ylab="Factor relevance", main="")

# Use n=2 PEER and save residues
peer.dir = paste0("data/gwas/rfu_normalize_peer2_", suffix)
# For matrixEQTL
peer_residuals = read.csv(paste(peer.dir, "/residuals.csv", sep=""), header=F)
rfu.normalized = read.table(paste0("data/gwas/subject_normalize_matrixeqtl_",suffix,".tsv"), 
                            header=T, sep='\t', row.names = 1)
colnames(peer_residuals) = colnames(rfu.normalized)
rownames(peer_residuals) = rownames(rfu.normalized)
write.table(peer_residuals, paste0("data/gwas/",suffix,"_residuals.tsv"), sep='\t', quote=F)


# For plots -------------------------
known.var = read_tsv("data/gwas/known_var.tsv")
ggplot(known.var, aes(x = log10(tcr.num))) +
  geom_histogram(fill = "#299D8F", color = "black") +
  labs(title = "", x = expression(log[10](TCR~count)), y = "Number of subjects") +
  theme_minimal()
ggsave("figs/tcrnum_samples.svg", width = 3.5, height = 3.5)

predicted.pop = read.table("data/qctmp/sel_rel_smps.txt", sep="\t", head=TRUE, comment.char = "")
predicted.pop[predicted.pop$Computed.population=="", "Computed.population"] = "Not assigned"
known.var.c = known.var %>% left_join(predicted.pop, by = c("IID" = "Subject"))

eigenval <- scan("./data/qctmp/sel_pca.eigenval")
pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval)*100)
plot(pve, ylab="Percent variance explained", xlab = "Principal component")

table(known.var$Computed.population)
table(known.var$sex) # female: 1, male: 2
subject.id.f = as_tibble(TEDDY.pheno) %>% filter(sample_id %in% known.var$IID) %>% pull(subject_id)
teddy.pheno.f = as_tibble(TEDDY.pheno) %>% filter(subject_id %in% subject.id.f, Assay.Type == "RNA-Seq")
mean(teddy.pheno.f$collinterval)/ 12
sd(teddy.pheno.f$collinterval) / 12

known.var.c = known.var %>% mutate(source = case_when(
  source == 1 ~ "Colorado",
  source == 2 ~ "Finland",
  source == 3 ~ "Georgia",
  source == 4 ~ "Germany",
  source == 5 ~ "Sweden",
  source == 6 ~ "Washington",
  TRUE ~ as.character(source)
))

# Plot for PCs and sources
ggplot(known.var.c, aes(PC1, PC2, color = source)) + geom_point() +
  coord_equal() + theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  labs(color="Subject source") + ggtitle("Genotype PC") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/rfuqtl_training_pc1_pc2.svg", width = 3.5, height = 3.5)

ggplot(known.var.c, aes(PC1, PC3, color = source)) + geom_point() +
  coord_equal() + theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) +
  labs(color="Subject source") + ggtitle("Genotype PC")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/rfuqtl_training_pc1_pc3.svg", width = 3.5, height = 3.5)

ggplot(known.var.c, aes(PC4, PC5, color = source)) + geom_point() +
  coord_equal() + theme_light() +
  xlab(paste0("PC4 (", signif(pve$pve[4], 3), "%)")) + ylab(paste0("PC5 (", signif(pve$pve[5], 3), "%)")) +
  labs(color="Subject source") + ggtitle("Genotype PC")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/rfuqtl_training_pc4_pc5.svg", width = 3.5, height = 3.5)

peer.dir = paste0("data/gwas/rfu_normalize_peer_", suffix)
peer_X = read.csv(paste(peer.dir, "/X.csv", sep=""), header=F)

peer_pca_cor = cor(t(peer_X[-1:-7,]), known.var %>% select(-IID), method = 'p', use="complete.obs")
rownames(peer_pca_cor) = as.numeric(rownames(peer_pca_cor)) - 7
peer_pca_cor2 = as.data.frame(peer_pca_cor) %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
ggplot(as.data.frame(peer_pca_cor2), aes(x = factor(rowname, levels=1:10), y = colname, fill = value)) + 
  geom_tile() + xlab("PEER factors") +
  ylab("Potential covariates") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                       midpoint = 0, 
                       guide = guide_colorbar(title = "Correlation")) +
  theme_minimal() +
  coord_equal()
ggsave("figs/rfuqtl_training_peer_cor.svg", width = 4.5, height = 4.5)

Alpha = read.csv(paste(peer.dir, "/Alpha.csv", sep=""), header=F)$V1

x_values <- 0:10
y_values <- 1.0 / Alpha[-1:-6]
data <- data.frame(x_values, y_values)
ggplot(data, aes(x=as_factor(x_values), y=y_values)) +
  geom_point() +
  labs(x = "Intercept and PEER Factors",
       y = "Factor relevance",
       title = "") +
  theme_minimal()
ggsave("figs/rfuqtl_training_peer_relevance.svg", width = 2.5, height = 2.5)
