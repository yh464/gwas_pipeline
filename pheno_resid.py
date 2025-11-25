#!/usr/bin/env python

def main(args):
    import os
    import pandas as pd
    import numpy as np
    from fnmatch import fnmatch
    from time import perf_counter as t
    from sklearn.linear_model import LinearRegression
    
    tic = t()
    def toc():
      print(f'Time = {t()-tic:.3f}')
    
    if not os.path.isdir(args.out): os.mkdir(args.out)
    
    # process covars
    dcov = pd.read_csv(args.cov, sep = '\s+')
    cov_out = dcov.iloc[:,:2]
    for i in range(2,dcov.shape[1]):
      tmp = dcov.iloc[:,i].values
      tmpname = dcov.columns.tolist()[i]
      tmpval = np.unique(tmp)
      for j in range(tmpval.size-1):
        cov_out.insert(cov_out.shape[1], column = f'{tmpname}_{tmpval[j]}',
                       value = (tmp==tmpval[j]).astype(np.int8))
    
    cov_out = cov_out.merge(pd.read_csv(args.qcov, sep = '\s+'), how = 'outer', on = ['FID', 'IID']).dropna()
    cov_out['const'] = np.ones(cov_out.shape[0])
    cov_collist = cov_out.columns.tolist()[2:]
    cov_id = cov_out[['FID','IID']]
    print('Covariates:')
    for x in cov_collist: print(x)
    
    toc()
    
    os.chdir(args._in)
    for x in os.listdir():
      if not fnmatch(x,'*.txt'): continue
      df = pd.read_csv(x, sep = '\s+').dropna()
      df = pd.merge(df, cov_id, on = ['FID','IID'], how = 'inner')
      out = df[['FID','IID']]
      ref = df[['FID','IID']]
      resid_list = []
      for y in df.columns.tolist()[2:]:
        endog = df[y].values
        exog = pd.merge(cov_out,ref, on = ['FID','IID'], how = 'inner').values
        md = LinearRegression().fit(exog, endog)
        resid = endog - md.predict(exog)
        resid_list.append(pd.Series(resid))
      out = pd.concat([out] + resid_list, axis = 1)
      out.columns = df.columns
      prefix = x.replace('.txt','')
      out.to_csv(f'{args.out}/{prefix}_resid.txt', index = False)
      print(f'Finished processing {x}')
      toc()
      del df, out, ref

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Computes correlational plots between IDPs and PGS (single file)')
    parser.add_argument('-i','--in', dest = '_in', help = 'Input directory',
      default = '../pheno/conn/')
    parser.add_argument('--cov',dest = 'cov', help = 'DISCRETE covariance file',
      default = '../params/discrete_covars.txt')
    parser.add_argument('--qcov',dest = 'qcov', help = 'QUANTITATIVE covariance file',
      default = '../params/quantitative_covars.txt')
    parser.add_argument('-o','--out', dest = 'out', help = 'Output directory',
      default = '../pheno/resid/')
    # always forces
    args=parser.parse_args()
    import os
    for arg in ['_in','out','cov','qcov']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    
    from ._utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in, __file__)
    proj.add_output(args._out, __file__)
    try: main(args)
    except: cmdhistory.errlog()