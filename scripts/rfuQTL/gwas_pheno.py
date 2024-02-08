#%%
import os
import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import scipy.stats as stats

os.chdir("/home/tany2/airr")
from scripts.gwas.normalize_function import prepare_expression

#%% gene expression data: G genes * N samples
fin = "results/rfu/subject_RFU.tsv"
fout = "data/qctmp/subject_rfu_matrixeqtl.tsv"

samples = pd.read_csv("data/qctmp/sel_rel.psam", sep='\t')

samples = samples[["#FID"]]

subject_RFU = pd.read_csv(fin, sep='\t')
subject_RFU.columns = subject_RFU.columns.str.split("_").str[2].astype(int)

teddy_pheno = pd.read_csv("data/rawdata/TEDDY.pheno.tsv", sep='\t')
wgs_subject = teddy_pheno.loc[teddy_pheno["sample_id"].str.startswith("WGS_"), ["sample_id", "subject_id"]]
wgs_subject.index = wgs_subject["subject_id"]
wgs_subject = wgs_subject["sample_id"]

subject_RFU.columns = np.vectorize(wgs_subject.get)(subject_RFU.columns.to_numpy())

subject_RFU = subject_RFU.T

subject_RFU["#FID"] = subject_RFU.index

wgs_rfu = pd.merge(samples, subject_RFU, "left")
print("missing value number", wgs_rfu.isna().sum().sum())

wgs_rfu = wgs_rfu.T

#%%
# Normalize
suffix = "tcrnum"
tpm_df = pd.read_csv("data/gwas/subject_rfu_matrixeqtl_"+suffix+".tsv", sep='\t', index_col=0)
print(f"{tpm_df.shape=}")

norm_df = prepare_expression(tpm_df)
print(f"{norm_df.shape=}")

# MatrixEQTL gene expression data: G genes * N samples
norm_df.to_csv("data/gwas/subject_normalize_matrixeqtl_"+suffix+".tsv", sep='\t')

# PEER experssion matrix: N samples * G genes
norm_df.T.to_csv("data/gwas/subject_normalize_peer_"+suffix+".tab", sep='\t', header=None, index=None)

