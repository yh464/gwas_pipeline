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

prefix = os.path.split(args._in)[1].replace('.txt','').split('/')[-1]
if not os.path.isfile(f'{args.out}/{prefix}_summary.csv') or args.force:
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
    
    # initialise output
    idp = np.tile(np.array(x_collist, dtype = 'U'), len(prs_collist)) # [idp1, idp2, ..., idp1, idp2, ...]
    prs = np.repeat(np.array(prs_collist, dtype = 'U'), len(x_collist)) # [prs1, ... prs1, prs2, ...]
    beta = []
    z = []
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
        beta.append(md.params[-1])
        z.append(md.tvalues[-1])
    
    summary = pd.DataFrame(dict(IDP = idp, PRS = prs, beta = beta, z = z)).dropna()
    summary['p'] = 1 - sts.chi2.cdf(summary['z']**2, df = 1)
    summary['q'] = sts.false_discovery_control(summary['p'])
    
    # summary table and tabular output
    tmp_beta = summary.pivot(index = 'IDP', columns = 'PRS', values = 'beta')
    tmp_z = summary.pivot(index = 'IDP', columns = 'PRS', values = 'z')
    tmp_beta.to_csv(f'{args.out}/{prefix}_beta.csv', sep = '\t',
                    index_label = False, index = True, header = True)
    tmp_z.to_csv(f'{args.out}/{prefix}_z.csv', sep = '\t', 
                 index_label = False, index = True, header = True)
    summary.to_csv(f'{args.out}/{prefix}_summary.csv', index = False, sep = '\t')

# plot figure
if not os.path.isfile(f'{args.out}/{prefix}.pdf') or args.force:
    import seaborn as sns
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    # Red-blue colour map
    cdict = dict(red = ((0,0,0),(1/2,1,1),(1,.8,.8)),
                 green = ((0,0,0),(1/2,1,1),(1,0,0)),
                 blue = ((0,.8,.8),(1/2,1,1),(1,0,0)))
    cmap_name = 'redblue'
    cmap = mpl.colors.LinearSegmentedColormap(cmap_name,cdict,1024)
    try:
      mpl.colormaps.register(cmap)
    except:
      pass
    
    summary = pd.read_table(f'{args.out}/{prefix}_summary.csv')
    # grid lines
    sns.set_theme(style = 'whitegrid')
    
    # point size
    summary['pt_size'] = (summary.q < 0.05).astype(float) + \
        (summary.p < 0.05).astype(float) + 1
        
    _, ax = plt.subplots(
        # figsize = (len(summary['PRS'].unique()) * 2.3,len(summary['IDP'].unique())))
        figsize = (7,len(summary['PRS'].unique())))
    
    beta_max = max([abs(summary['beta'].max()), abs(summary['beta'].min())])
    sns.scatterplot(
        summary,
        x = 'PRS', y = 'IDP',
        hue = 'beta', palette = 'redblue', hue_norm = (-beta_max, beta_max),
        size = 'pt_size', sizes = (50, 400),
        edgecolor = '.7',
        legend = False, ax = ax
        )
    # sns.scatterplot(
    #     summary,
    #     x = 'PRS', y = 'IDP',
    #     hue = 'z', palette = 'redblue', hue_norm = (-2.5,2.5),
    #     size = 'pt_size', sizes = (50, 400),
    #     edgecolor = '.7',
    #     legend = False, ax = ax[1]
    #     )
    
    for _, spine in ax.spines.items():
          spine.set_visible(False)
    for label in ax.get_xticklabels():
          label.set_rotation(90)
    plt.savefig(f'{args.out}/{prefix}.png')
    plt.savefig(f'{args.out}/{prefix}.pdf')
    plt.close()
toc = time.perf_counter()-tic
print(f'Finished processing {prefix}. time = {toc:.3f} seconds')