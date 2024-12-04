#!/usr/bin/env python

def main(args):
    import pandas as pd
    from fnmatch import fnmatch
    import numpy as np
    
    f = open(args._in).read().splitlines()
    pf = []
    for x in f:
      pf.append(x.split('/')[-1].split('.')[0])
    
    df = pd.read_csv(f[0], sep = '\s+')
    for i in range(1,len(f)):
      df = df.merge(pd.read_csv(f[i], sep = '\s+'), on = ['FID','IID'],
                    suffixes = (pf[i-1],pf[i]))
    
    for i in df.columns:
      if fnmatch(i, '*frac*') or fnmatch(i, '*abs*'):
        df.drop(i, axis = 'columns', inplace = True)
        continue
      if sum(df[i].isna()) > 0.5 * len(df[i]): 
        df.drop(i, axis = 'columns', inplace = True)
    c = df.iloc[:,2:].corr() # remove columns FID and IID, the corr() fcn ignores NaN
    c = np.abs(c)
    l,_ = np.linalg.eig(c)
    n = l.size
    neff_abs = n * (1 - (n-1) * np.var(l) / n**2)
    l[l<0] = 0
    neff_pos = n * (1 - (n-1) * np.var(l) / n**2)
    
    with open(args.out,'w') as f:
      print(f'Effective # variables by old method Nyholt DR (2004): {neff_abs}', file = f)
      print(f'Bonferroni corrected threshold: {1-0.95**(1/neff_abs)}', file = f)
      print(file = f)
      print(f'Effective # variables by new method Nyholt DR (2004): {neff_pos}', file = f)
      print(f'Bonferroni corrected threshold: {1-0.95**(1/neff_pos)}', file = f)
      
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
      description = 'This programme calculates the effective number of indep. variables.')
    parser.add_argument('-i','--in', dest = '_in', help = 'List containing all phenotypes',
      default = '../params/pheno_list.txt')
    parser.add_argument('-o', '--out', dest = 'out', help = 'output file',
      default = '../pheno/neff.txt')
    # always overwrites
    args = parser.parse_args()
    import os
    for arg in ['_in','out']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in, __file__)
    proj.add_output(args.out, __file__)
    try: main(args)
    except: cmdhistory.errlog()