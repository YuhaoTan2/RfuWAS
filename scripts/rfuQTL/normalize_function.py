import os
import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import scipy.stats as stats

##### RFU normalization

def normalize_quantiles(df):
    """
    Quantile normalization to the average empirical distribution
    Note: replicates behavior of R function normalize.quantiles
          from library("preprocessCore")
    Reference:
     [1] Bolstad et al., Bioinformatics 19(2), pp. 185-193, 2003
    Adapted from https://github.com/andrewdyates/quantile_normalize
    """
    M = df.values.copy()

    Q = M.argsort(axis=0)
    m,n = M.shape

    # compute quantile vector
    quantiles = np.zeros(m)
    for i in range(n):
        quantiles += M[Q[:,i],i]
    quantiles = quantiles / n

    for i in range(n):
        # Get equivalence classes; unique values == 0
        dupes = np.zeros(m, dtype=int)
        for j in range(m-1):
            if M[Q[j,i],i] == M[Q[j+1,i],i]:
                dupes[j+1] = dupes[j]+1

        # Replace column with quantile ranks
        M[Q[:,i],i] = quantiles

        # Average together equivalence classes
        j = m-1
        while j >= 0:
            if dupes[j] == 0:
                j -= 1
            else:
                idxs = Q[j-dupes[j]:j+1,i]
                M[idxs,i] = np.median(M[idxs,i])
                j -= 1 + dupes[j]
        assert j == -1

    return pd.DataFrame(M, index=df.index, columns=df.columns)


def inverse_normal_transform(M):
    """Transform rows to a standard normal distribution"""
    if isinstance(M, pd.Series):
        r = stats.rankdata(M)
        return pd.Series(stats.norm.ppf(r/(M.shape[0]+1)), index=M.index, name=M.name)
    else:
        R = stats.rankdata(M, axis=1)  # ties are averaged
        Q = stats.norm.ppf(R/(M.shape[1]+1))
        if isinstance(M, pd.DataFrame):
            Q = pd.DataFrame(Q, index=M.index, columns=M.columns)
        return Q


def prepare_expression(tpm_df, sample_frac_threshold=0.2,
                       count_threshold=6, tpm_threshold=0, mode='qn'):
    """
    Genes are filtered using the following expression thresholds:
      TPM >= tpm_threshold in >= sample_frac_threshold * samples
      read counts >= count_threshold in sample_frac_threshold * samples

    The filtered counts matrix is then normalized using:
      TMM (mode='tmm'; default) or
      quantile normalization (mode='qn')
    Previous default: count_threshold=6, tpm_threshold=0.1
    """

    # expression thresholds
    ns = tpm_df.shape[1]
    mask = (
        (np.sum(tpm_df > tpm_threshold, axis=1) >= sample_frac_threshold * ns)
        ).values
    print("Filter out genes:", np.where(mask==0)[0], np.where(mask==0)[0].shape)
    # apply normalization
    if mode.lower() == 'tmm':
        raise NotImplementedError('TMM normalization not implemented')
    elif mode.lower() == 'qn':
        #qn_df =  normalize_quantiles(tpm_df)
        qn_df =  normalize_quantiles(tpm_df.loc[mask])
        norm_df =  inverse_normal_transform(qn_df)
    else:
        raise ValueError(f'Unsupported mode {mode}')

    return norm_df
