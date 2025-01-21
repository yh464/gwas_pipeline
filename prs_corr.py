#!/usr/bin/env python3
'''
Computes correlational plots between IDPs and PGS
'''

import argparse
parser = argparse.ArgumentParser(description='Computes correlational plots between IDPs and PGS (single file)')
parser.add_argument('-i','--in', dest = '_in', help = 'Input file')
parser.add_argument('-p','--prs', dest = 'prs', help = 'PRS score directory',
  default = '../prs/prs_score/')
parser.add_argument('--cov',dest = 'cov', help = 'DISCRETE covariance file',
  default = '../params/discrete_covars.txt')
parser.add_argument('--qcov',dest = 'qcov', help = 'QUANTITATIVE covariance file',
  default = '../params/quantitative_covars.txt')
parser.add_argument('-o','--out', dest = 'out', help = 'Output directory',
  default = '../prs/prs_corr/')
parser.add_argument('-f','--force', dest = 'force', help = 'force overwrite',
                    action = 'store_true', default = False)
args=parser.parse_args()
import os
for arg in ['_in','out','cov','qcov','prs']:
    exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
import time
tic = time.perf_counter()

from fnmatch import fnmatch
import pandas as pd

prefix = os.path.basename(args._in).replace('.txt','')
if not os.path.isfile(f'{args.out}/{prefix}_summary.txt') or args.force:
    import numpy as np
    from statsmodels.api import OLS
    import scipy.stats as sts
    
    # process covars
    dcov = pd.read_csv(args.cov, sep = '\s+')
    cov_out = dcov.iloc[:,:2]
    for i in range(2,dcov.shape[1]):
      tmp = dcov.iloc[:,i].values
      tmpname = dcov.columns.tolist()[i]
      tmpval = np.unique(tmp)
      for j in range(tmpval.size-1):
        cov_out[f'{tmpname}_{tmpval[j]}'] = (tmp==tmpval[j]).astype(np.int8)
    
    cov_out = cov_out.merge(pd.read_csv(args.qcov, sep = '\s+'), how = 'outer', on = ['FID', 'IID']).dropna()
    cov_out['const'] = np.ones(cov_out.shape[0])
    cov_collist = cov_out.columns.tolist()[2:]
    
    # process PRS
    os.chdir(args.prs)
    prslist = []
    prs_collist = []
    for x in os.listdir():
      if not fnmatch(x, '*.txt'): continue
      prs = x.replace('.txt','')
      prs_collist.append(prs)
      tmp = pd.read_csv(x,sep = '\s+').iloc[:,[0,1,3]]
      tmp.columns = ['FID','IID',prs]
      prslist.append(tmp)
    
    prsdf = prslist[0]
    for x in range(1,len(prslist)):
      prsdf = prsdf.merge(prslist[x], how = 'outer', on = ['FID','IID'])
    toc = time.perf_counter() - tic
    print(f'Finished processing PRS input. time = {toc:.3f} seconds')
    
    # read input
    tmp = pd.read_csv(args._in,sep = '\s+')
    x_list = tmp.columns.tolist()[2:]
    
    x_collist = []
    for xtmp in x_list:
      if fnmatch(xtmp,'*fra*') or fnmatch(xtmp, '*abs*') or fnmatch(xtmp,'bet*'):
        continue
      x_collist.append(xtmp)
    
    tmp = tmp.merge(cov_out, how = 'inner', on = ['FID', 'IID'])
    tmp = tmp.merge(prsdf, how = 'inner', on = ['FID', 'IID'])
    
    summary = []
    for i in range(len(prs_collist)):
      tmp_collist = cov_collist
      tmp_collist.append(prs_collist[i])
      
      # linear model using cov, const and prs as exog variable
      for j in range(len(x_collist)):
        endog = tmp[x_collist[j]].values
        endog = (endog - endog.mean())/endog.std()
        # normalising to z score might be misleading - consider commenting out
        
        exog = tmp[tmp_collist].values
        md = OLS(endog, exog).fit()
        summary.append(
            pd.DataFrame(dict(
                group1 = prefix,
                pheno1 = x_collist[j],
                group2 = ['PRS'],
                pheno2 = prs_collist[i],
                beta = md.params[-1],
                z = md.tvalues[-1]
                ))
            )
        
    summary = pd.concat(summary).dropna()
    summary['p'] = 1 - sts.chi2.cdf(summary['z']**2, df = 1)
    summary['q'] = sts.false_discovery_control(summary['p'])
    
    # summary table and tabular output
    tmp_beta = summary.pivot(index = 'pheno1', columns = 'pheno2', values = 'beta')
    tmp_beta.columns.name = None; tmp_beta.index.name = None
    tmp_beta.to_csv(f'{args.out}/{prefix}_beta.txt', sep = '\t',
                    index_label = False, index = True, header = True)
    summary.to_csv(f'{args.out}/{prefix}_summary.txt', index = False, sep = '\t')