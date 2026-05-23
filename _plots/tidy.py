#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2026-05-23

Data tidying for plotting
'''
import pandas as pd
import numpy as np
from scipy.stats import false_discovery_control as fdr

def capitalise(series):
    out = []
    for x in range(len(series)):
        tmp = series.iloc[x]
        if type(tmp) != str: out.append(tmp)
        elif len(tmp) < 1: out.append(tmp)
        elif all([x.isupper() or not x.isalpha() for x in tmp]): out.append(tmp) # if all uppercase, keep as is
        else:
            tmp = tmp[0].upper() + tmp[1:]
            out.append(tmp)
    return out

def get_fdr(df, group_by = [0,1], sig_col = None, p_threshold: list[float] = []):
    if 'fdr' in df.columns: df['q'] = df['fdr']
    if 'FDR' in df.columns: df['q'] = df['fdr']
    if 'P' in df.columns: df['p'] = df['P']
    
    if sig_col is not None and sig_col in df.columns:
        sig_col_copy = df[sig_col].copy()
        df = df.assign(Significance = 'not significant')
        df.loc[sig_col_copy == True, 'Significance'] = 'significant'
        sig_label = 'significant'
    
    elif 'p' in df.columns and len(p_threshold) == 0:
        if not 'q' in df.columns:
            df = df.assign(q = np.nan)
            if len(group_by) == 0: # global FDR
                df.loc[~np.isnan(df['p']),'q'] = fdr(df.loc[~np.isnan(df['p']),'p'])
            else:
                for _, df_group in df.groupby(group_by):
                    df.loc[df.index.isin(df_group.index) & ~np.isnan(df_group['p']),'q'] = \
                        fdr(df_group.loc[~np.isnan(df_group['p']),'p'])
        df = df.assign(Significance = 'not significant')
        df.loc[df['p'] < 0.05, 'Significance'] = 'nominal'
        df.loc[df['q'] < 0.05, 'Significance'] = 'FDR-corrected'
        sig_label = 'FDR-corrected'

    elif 'p' in df.columns and len(p_threshold) > 0:
        p_threshold = sorted(p_threshold, reverse = True)
        df = df.assign(Significance = 'not significant')
        for p_thr in p_threshold:
            df.loc[df['p'] < p_thr, 'Significance'] = f'p < {p_thr}'
        sig_label = f'p < {p_threshold[-1]}'

    else:
        df = df.assign(Significance = 'NA')
        sig_label = ''

    return df, sig_label