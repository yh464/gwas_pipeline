#!/usr/bin/env python3
'''
Computes correlational plots between IDPs and PGS
'''

import argparse
parser = argparse.ArgumentParser(description='Computes correlational plots between IDPs and PGS (single file)')
parser.add_argument('-i','--in', dest = '_in', help = 'Input file')
parser.add_argument('-p','--prs', dest = 'prs', help = 'PRS score directory',
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
for arg in ['_in','out','dcov','qcov','prs']:
    exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
import time
tic = time.perf_counter()

# from fnmatch import fnmatch
import pandas as pd

prefix = os.path.basename(args._in).replace('.txt','')
if not os.path.isfile(f'{args.out}/{prefix}_summary.txt') or args.force:
    import numpy as np
    from statsmodels.api import OLS
    import scipy.stats as sts
    
    # process covars
    dcov = pd.read_table(args.dcov, index_col = ['FID','IID'], sep = '\\s+').astype(str)
    cov_out = pd.get_dummies(dcov, prefix = dcov.columns, drop_first = True)
    # for i in range(dcov.shape[1]):
    #   tmp = dcov.iloc[:,i].values
    #   tmpname = dcov.columns.tolist()[i]
    #   tmpval = np.unique(tmp)
    #   for j in range(tmpval.size-1):
    #     cov_out[f'{tmpname}_{tmpval[j]}'] = (tmp==tmpval[j]).astype(np.int8)
    
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
    os.chdir(args.prs)
    prs_list = [x.replace('.txt','') for x in os.listdir()]
    prs_list = np.array(prs_list, dtype = 'U')
    prs_list = np.unique(prs_list)
    for prs in prs_list:
        if prs == 'qc': continue
        if not os.path.isdir(prs) and not os.path.isfile(f'{prs}.txt'): continue
        print(prs)
        # sum PRS across all chromosomes
        if (not os.path.isfile(f'{args.prs}/{prs}.txt') or \
            args.force) and os.path.isdir(f'{args.prs}/{prs}'):
            chrom = []
            for j in range(22):
                df = pd.read_table(f'{prs}/{prs}.chr{j+1}.sscore').sort_values('IID').drop_duplicates()
                col = df.columns.tolist()
                col[0] = 'FID'
                df.columns = col
                df = df.set_index(['FID','IID'])
                chrom.append(df.iloc[:,-1].to_frame())
            chrom = pd.concat(chrom, axis = 1)
            out = pd.DataFrame(index = chrom.index, columns = [])
            out['score_total'] = chrom.sum(axis = 1)
            out['score_norm'] = out.score_total/out.score_total.std()
            out.to_csv(f'{args.prs}/{prs}.txt', index = True, sep = '\t')
        else:
            out = pd.read_table(f'{args.prs}/{prs}.txt', index_col = ['FID','IID'])
        print(f'    {args.prs}/{prs}.txt')
        
        tmp = out[['score_norm']]
        tmp.columns = [prs]
        tmp_merge = pd.concat([phen_cov.copy(),tmp], axis = 1, join = 'inner')
        print(f'    Sample size: {tmp_merge.shape[0]}')
        for phen in phen_list:
            tmp = tmp_merge[cov_list + [prs] + [phen]].dropna().astype(float)
            endog = tmp[phen].values
            exog = tmp[cov_list + [prs]].values
            try:
                md = OLS(endog, exog).fit()
                summary.append(
                    pd.DataFrame(dict(
                        group1 = prefix,
                        pheno1 = phen,
                        group2 = ['PRS'],
                        pheno2 = prs,
                        beta = md.params[-1],
                        z = md.tvalues[-1],
                        n = tmp.shape[0]                        
                        ))
                    )
            except:
                print(f'        {phen} correlation failed with {prs}')
    summary = pd.concat(summary).dropna()
    summary['p'] = 1 - sts.chi2.cdf(summary['z']**2, df = 1)
    summary['q'] = sts.false_discovery_control(summary['p'])
    
    # summary table and tabular output
    tmp_beta = summary.pivot(index = 'pheno1', columns = 'pheno2', values = 'beta')
    tmp_beta.columns.name = None; tmp_beta.index.name = None
    tmp_beta.to_csv(f'{args.out}/{prefix}_beta.txt', sep = '\t',
                    index_label = False, index = True, header = True)
    summary.to_csv(f'{args.out}/{prefix}_summary.txt', index = False, sep = '\t')