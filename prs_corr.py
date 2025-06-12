#!/usr/bin/env python3
'''
Computes correlational plots between IDPs and PGS
'''

import argparse
parser = argparse.ArgumentParser(description='Computes correlational plots between IDPs and PGS (batch)')
parser.add_argument('-i','--in', dest = '_in', help = 'Input directory',
  default = '../pheno/ukb/')
parser.add_argument('-p','--prs', dest = 'prs', help = 'PRS to select', nargs = '*',
  default = ['disorders','disorders_subtypes'])
parser.add_argument('--prsdir', dest = 'prsdir', help = 'PRS score directory',
  default = '../prs/prs_score/')
parser.add_argument('--dcov',dest = 'dcov', help = 'DISCRETE covariance file',
  default = '../params/discrete_covars.txt')
parser.add_argument('--qcov',dest = 'qcov', help = 'QUANTITATIVE covariance file',
  default = '../params/quantitative_covars.txt')
parser.add_argument('-o','--out', dest = 'out', help = 'Output directory',
  default = '../prs/prs_corr/')
parser.add_argument('-f','--force', dest = 'force', help = 'force overwrite',
                    action = 'store_true', default = False)
args=parser.parse_args()
import os
for arg in ['_in','out','dcov','qcov','prsdir']:
    setattr(args, arg, os.path.realpath(getattr(args, arg)))

import time
tic = time.perf_counter()
import pandas as pd

prefix = os.path.basename(args._in).replace('.txt','')
if not os.path.isfile(f'{args.out}/{prefix}_summary.txt') or args.force:
    import pingouin as pg
    import numpy as np
    from statsmodels.api import OLS
    import scipy.stats as sts
    
    # process covars
    dcov = pd.read_table(args.dcov, index_col = ['FID','IID'], sep = '\\s+').astype(str)
    cov_out = pd.get_dummies(dcov, prefix = dcov.columns, drop_first = True).astype(int)
    cov_out = pd.concat([cov_out,pd.read_table(args.qcov, index_col = ['FID','IID'], sep = '\\s+')], 
                        axis = 1)
    cov_list = cov_out.columns.tolist()
    
    # read input
    phen = pd.read_table(args._in, index_col = ['FID','IID'], sep = '\\s+')
    # standardise phenotypes
    print('Phenotypes:')
    phen_list = []
    for x in phen.columns: 
        print(f'    {x}')
        try:
            phen.loc[~phen[x].isna(), x] -= phen.loc[~phen[x].isna(), x].mean()
            phen.loc[~phen[x].isna(), x] /= phen.loc[~phen[x].isna(), x].std()
            phen_list.append(x)
        except:
            print(f'    {x} is not quantitative, dropping')
            phen.drop(x, inplace = True, axis = 1)
    phen_cov = pd.concat([phen,cov_out], axis = 1)
    print(f'Total sample size: {phen_cov.shape[0]}\n')
    
    summary = []
    # correlate for each PRS
    for p in args.prs:
        prs_list = []
        for x in os.listdir(f'{args.prsdir}/{p}'):
            if x == 'qc': continue
            if os.path.isfile(f'{args.prsdir}/{p}/{x}.txt') or os.path.isdir(f'{args.prsdir}/{p}/{x}'):
                prs_list.append(x)
        
        for prs in prs_list:
            print(prs)
            # sum PRS across all chromosomes
            if (not os.path.isfile(f'{args.prsdir}/{p}/{prs}.txt') or \
                args.force) and os.path.isdir(f'{args.prsdir}/{p}/{prs}'):
                chrom = []
                for j in range(22):
                    df = pd.read_table(f'{args.prsdir}/{p}/{prs}/{prs}.chr{j+1}.sscore'
                        ).sort_values('IID').drop_duplicates()
                    col = df.columns.tolist()
                    col[0] = 'FID'
                    df.columns = col
                    df = df.set_index(['FID','IID'])
                    chrom.append(df.iloc[:,-1].to_frame())
                chrom = pd.concat(chrom, axis = 1)
                out = pd.DataFrame(index = chrom.index, columns = [])
                out['score_total'] = chrom.sum(axis = 1)
                out['score_norm'] = out.score_total/out.score_total.std()
                out.to_csv(f'{args.prsdir}/{p}/{prs}.txt', index = True, sep = '\t')
            else:
                out = pd.read_table(f'{args.prsdir}/{p}/{prs}.txt', index_col = ['FID','IID'])
            print(f'    {args.prsdir}/{p}/{prs}.txt')
            
            tmp = out[['score_norm']]
            tmp.columns = ['prs']
            tmp_merge = pd.concat([phen_cov.copy(),tmp], axis = 1, join = 'inner')
            print(f'    Sample size: {tmp_merge.shape[0]}')
            for phen in phen_list:
                try:
                    tmp = tmp_merge[cov_list + ['prs',phen]].dropna()
                    partcorr = pg.partial_corr(tmp, x = 'prs', y = phen, y_covar = cov_list)
                    endog = tmp[phen].values; exog = tmp[cov_list + ['prs']].values
                    md = OLS(endog, exog).fit()
                    ols_beta = md.params[-1]; ols_p = 1-sts.chi2.cdf(md.tvalues[-1]**2, df=1)
                    summary.append(
                        pd.DataFrame(dict(
                            group1 = prefix,
                            pheno1 = phen,
                            group2 = [p],
                            pheno2 = prs,
                            r = partcorr['r'].values[0],
                            p = partcorr['p-val'].values[0],
                            ols_beta = ols_beta, ols_p = ols_p,
                            n = partcorr['n'].values[0]                       
                            ))
                        )
                except:
                    print(f'        {phen} correlation failed with {prs}')
    summary = pd.concat(summary).dropna()
    summary['q'] = sts.false_discovery_control(summary['p'])
    summary['ols_q'] = sts.false_discovery_control(summary['ols_p'])
    
    # summary table and tabular output
    tmp_beta = summary.pivot(index = 'pheno1', columns = 'pheno2', values = 'r')
    tmp_beta.columns.name = None; tmp_beta.index.name = None
    tmp_beta.to_csv(f'{args.out}/{prefix}_beta.txt', sep = '\t',
                    index_label = False, index = True, header = True)
    summary.to_csv(f'{args.out}/{prefix}_summary.txt', index = False, sep = '\t')