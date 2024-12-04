#!/usr/bin/env python3
def main(args):
  import pandas as pd
  import numpy as np
  from statsmodels.api import OLS, Logit
  # from statsmodels.genmod.families.family import Binomial
  from fnmatch import fnmatch
  from time import perf_counter as t

  tic = t()
  
  import os
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
  
  # process input file
  fin = args._in
  prefix = os.path.basename(fin).replace('.txt','')
  # dfin = pd.read_csv(fin) # if reading from pheno/resid
  dfin = pd.read_csv(fin, sep = '\s+') # if reading from pheno/conn
  x_collist = dfin.columns.tolist()[2:]
  print(x_collist)
  
  # progress check
  bin_done = (os.path.isfile(f'{args.out}/{prefix}_ukbbinary_beta.csv') and 
              os.path.isfile(f'{args.out}/{prefix}_ukbbinary_z.csv'))
  quant_done = (os.path.isfile(f'{args.out}/{prefix}_ukbquant_beta.csv') and 
              os.path.isfile(f'{args.out}/{prefix}_ukbquant_z.csv'))
  
  if not (bin_done and quant_done) or args.force:
    # process covars
    dcov = pd.read_csv(args.cov, sep = '\s+')
    cov = dcov.iloc[:,:2]
    for i in range(2,dcov.shape[1]):
      tmp = dcov.iloc[:,i].values
      tmpname = dcov.columns.tolist()[i]
      tmpval = np.unique(tmp)
      for j in range(tmpval.size-1):
        cov[f'{tmpname}_{tmpval[j]}'] = (tmp==tmpval[j]).astype(np.int8)
  
    cov = cov.merge(pd.read_csv(args.qcov, sep = '\s+'), how = 'outer', on = ['FID', 'IID']).dropna()
    cov['const'] = np.ones(cov.shape[0])
    cov_collist = cov.columns.tolist()[2:]
    print(cov_collist)
    
    # process UKB phenotypes, binary and quantiative data processed separately (!)
    os.chdir(args.ukb)
    binlist = []
    quantlist = []
    for x in os.listdir():
      if fnmatch(x, '*binary.txt'):
        tmp = pd.read_csv(x,sep = '\s+')
        binlist.append(tmp)
      if fnmatch(x, '*quant.txt'):
        tmp = pd.read_csv(x,sep = '\s+')
        quantlist.append(tmp)
  
  if not bin_done or args.force:
    dfbin = binlist[0]
    for x in range(1,len(binlist)):
      dfbin = dfbin.merge(binlist[x], how = 'outer', on = ['FID','IID'])
    bin_collist = dfbin.columns.tolist()[2:]
    print(bin_collist)
  
  if not quant_done or args.force:
    dfquant = quantlist[0]
    for x in range(1,len(quantlist)):
      dfquant = dfquant.merge(quantlist[x], how = 'outer', on = ['FID','IID'])  
    quant_collist = dfquant.columns.tolist()[2:]
    print(quant_collist)
  toc = t() - tic
  print(f'Finished processing input. time = {toc:.3f} seconds')
  
  # correlations
  if not bin_done or args.force:
    master_dfbin = dfin.merge(dfbin, how = 'inner', on = ['FID', 'IID']).merge(
      cov, how = 'inner', on = ['FID', 'IID'])
    bin_beta = np.zeros((len(x_collist),len(bin_collist)))
    bin_z = bin_beta.copy()
    
    for i in range(len(x_collist)):
      for j in range(len(bin_collist)):
        tmp_collist = cov_collist.copy()
        tmp_collist.append(x_collist[i])
        tmp_collist.append(bin_collist[j])
        print(f'[{x_collist[i]}, {bin_collist[j]}]')
        tmp = master_dfbin[tmp_collist].dropna() # important as UKB data is not full!
        
        # endog = tmp[bin_collist[j]]
        # exog = tmp[tmp_collist[:-1]]
        
        endog = tmp[x_collist[i]]
        tmp_collist.remove(x_collist[i])
        exog = tmp[tmp_collist]
        try:
          md = Logit(endog, exog).fit(disp=0)
          bin_beta[i,j] = md.params[-1]
          bin_z[i,j] = md.tvalues[-1]
        except:
          bin_beta[i,j] = np.nan
          bin_z[i,j] = np.nan
        print(f'beta = {bin_beta[i,j]}, z = {bin_z[i,j]}')
    
    # output for bin
    bin_beta = pd.DataFrame(bin_beta, index = x_collist, columns = bin_collist)
    bin_beta.to_csv(f'{args.out}/{prefix}_ukbbinary_beta.csv', sep = '\t', index = True)
    bin_z = pd.DataFrame(bin_z, index = x_collist, columns = bin_collist)
    bin_z.to_csv(f'{args.out}/{prefix}_ukbbinary_z.csv', sep = '\t', index = True)
  else:
    bin_beta = pd.read_table(f'{args.out}/{prefix}_ukbbinary_beta.csv', index_col = 0)
    bin_z = pd.read_table(f'{args.out}/{prefix}_ukbbinary_z.csv', index_col = 0)
  toc = t() - tic
  print(f'Finished binary correlations. time = {toc:.3f} seconds')
  
  if not quant_done or args.force:
    master_dfquant = dfin.merge(dfquant, how = 'inner', on = ['FID', 'IID']).merge(
      cov, how = 'inner', on = ['FID', 'IID'])
    quant_beta = np.zeros((len(x_collist),len(quant_collist)))
    quant_z = quant_beta.copy()
    
    for i in range(len(x_collist)):
      for j in range(len(quant_collist)):
        tmp_collist = cov_collist.copy()
        tmp_collist.append(x_collist[i])
        tmp_collist.append(quant_collist[j])
        print(f'[{x_collist[i]}, {quant_collist[j]}]')      
        tmp = master_dfquant[tmp_collist].dropna() # important as UKB data is not full!
        
        # endog = tmp[quant_collist[j]]
        # exog = tmp[tmp_collist[:-1]]
        
        endog = tmp[x_collist[i]]
        tmp_collist.remove(x_collist[i])
        exog = tmp[tmp_collist]
        
        md = OLS(endog, exog).fit()
        quant_beta[i,j] = md.params[-1]
        quant_z[i,j] = md.tvalues[-1]
        print(f'beta = {quant_beta[i,j]}, z = {quant_z[i,j]}')
    tmp_collist = tmp_collist[:-1]
    toc = t() - tic
    
    # output
    quant_beta = pd.DataFrame(quant_beta, index = x_collist, columns = quant_collist)
    quant_beta.to_csv(f'{args.out}/{prefix}_ukbquant_beta.csv', sep = '\t', index = True)
    quant_z = pd.DataFrame(quant_z, index = x_collist, columns = quant_collist)
    quant_z.to_csv(f'{args.out}/{prefix}_ukbquant_z.csv', sep = '\t', index = True)
  else:
    quant_beta = pd.read_table(f'{args.out}/{prefix}_ukbquant_beta.csv', index_col = 0)
    quant_z = pd.read_table(f'{args.out}/{prefix}_ukbquant_z.csv', index_col = 0)
  print(f'Finished quantitative variable correlations. time = {toc:.3f} seconds')
  
  if not os.path.isfile(f'{args.out}/{prefix}_ukbbinary.png') or args.force:
    _, ax = plt.subplots(1,2,figsize = (bin_beta.shape[1]*2.3, bin_beta.shape[0]*.8))  
    sns.heatmap(bin_beta, annot = True, fmt = '.4f', ax = ax[0], cmap = 'redblue')
    sns.heatmap(bin_z, annot = True, fmt = '.4f',
                ax = ax[1], cmap = 'redblue', vmin = -2.5, vmax = 2.5)
    plt.savefig(f'{args.out}/{prefix}_ukbbinary.png')
    plt.close()
    
  if not os.path.isfile(f'{args.out}/{prefix}_ukbquant.png') or args.force:
    _, ax = plt.subplots(1,2,figsize = (quant_beta.shape[1]*2.3, quant_beta.shape[0]*.8))  
    sns.heatmap(quant_beta, annot = True, fmt = '.4f', ax = ax[0], cmap = 'redblue')
    sns.heatmap(quant_z, annot = True, fmt = '.4f',
                ax = ax[1], cmap = 'redblue', vmin = -2.5, vmax = 2.5)
    plt.savefig(f'{args.out}/{prefix}_ukbquant.png')
    plt.close()
  
if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Computes correlational plots between IDPs and PGS (single file)')
  parser.add_argument('-i','--in', dest = '_in', help = 'Input IDP file')
  parser.add_argument('-u','--ukb', dest = 'ukb', help = 'UKB phenotype directory',
    default = '../pheno/ukb/')
  parser.add_argument('--cov',dest = 'cov', help = 'DISCRETE covariance file',
    default = '../params/discrete_covars.txt')
  parser.add_argument('--qcov',dest = 'qcov', help = 'QUANTITATIVE covariance file',
    default = '../params/quantitative_covars.txt')
  parser.add_argument('-o','--out', dest = 'out', help = 'Output directory',
    default = '../gene_corr/pcorr/')
  parser.add_argument('-f','--force', dest = 'force', help = 'force overwrite',
                      action = 'store_true', default = False)
  args=parser.parse_args()
  
  import os
  for arg in ['_in','out','ukb','cov','qcov']:
      exec(f'args.{arg} = os.path.realpath(args.{arg})')
  main(args)