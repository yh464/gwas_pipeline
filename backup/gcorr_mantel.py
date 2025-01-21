#!/usr/bin/env python3
'''
This script performs the Mantel test for genetic / phenotypic similarity matrices
'''

def main(args):
    import numpy as np
    from scipy.stats import pearsonr
    import pandas as pd
    import time
    
    rng = np.random.Generator(np.random.MT19937(seed = 464))
    
    def mantelr(m1, m2):
      n = m1.shape[0]
      if m1.shape[1] != n or (m2.shape[0] != n or m2.shape[1] != n):
        raise ValueError('Must compare two SQUARE matrices with the same size')
      lt = np.tri(n, k=-1).astype('?')
      res = pearsonr(m1[lt], m2[lt])
      return res.statistic
      
    for x in args._in:
      df = pd.read_csv(x, index_col = 0)
      n = df.shape[0]
      if df.shape[1] != n:
        print(f'WARNING: {x} is not symmetric, skipping')
        continue
      
      lt = np.tri(n).astype('?')
      
      pcorr = df.values.copy()
      gcorr = df.values.copy()
      
      pcorr[lt.T] = 0; pcorr += pcorr.T                                            # upper triangle = 0, take lower triangle
      gcorr[lt] = 0; gcorr += gcorr.T
      
      
      mantel_r = mantelr(pcorr, gcorr)
      print(mantel_r)
      
      # permutation test, permute rows/columns
      null_dist = np.zeros(args.nperm)
      
      tic = time.perf_counter()
      for i in range(args.nperm):
        perm = rng.permutation(n)
        null_dist[i] = mantelr(pcorr, gcorr[perm,:][:,perm])
        if i % 5000 == 0:
          toc = time.perf_counter() - tic
          print(f'{toc:.2f} seconds, {i}/{args.nperm}')
      
      # parse p value
      p = (null_dist > abs(mantel_r)).sum() + (null_dist < -abs(mantel_r)).sum()
      p /= args.nperm
      
      # output log
      log = x.replace('.csv', '.mantel.log')
      f = open(log, 'w')
      print(f'Mantel R = {mantel_r:.5f}', file = f)
      print(f'p value = {p:.4e}', file = f)
      for _ in range(5): print(file = f)
      print('Mantel permutation cache:', file = f)
      for i in null_dist: print(f'{i:.4f}', file = f)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
      description='Mantel test script for genetic/phenotypic correlations')
    parser.add_argument('-i','--in', dest = '_in', help = 'input files', nargs = '*')
    parser.add_argument('--nperm', dest = 'nperm', help = '# permutations', type = int, 
                        default = 10000)
    args = parser.parse_args()
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    for i in args._in: proj.add_input(i,__file__)
    try: main(args)
    except: cmdhistory.errlog()