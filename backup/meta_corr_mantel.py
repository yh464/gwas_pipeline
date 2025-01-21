#!/usr/bin/env python3
    
def main(args):
    import numpy as np
    from scipy.stats import pearsonr
    import pandas as pd
    import time
    
    rng = np.random.Generator(np.random.MT19937(seed = 464))
    
    def mantelr(m1,m2):
      if m1.shape != m2.shape: raise ValueError('Two matrices need to be the same shape')
      return pearsonr(m1.reshape(m1.size),m2.reshape(m2.size)).statistic
    
    pf1 = args._in[0].split('/')[-1].replace('.csv','')
    pf2 = args._in[1].split('/')[-1].replace('.csv','')
    
    m1 = pd.read_csv(args._in[0], index_col = 0).values()
    m2 = pd.read_csv(args._in[1], index_col = 0).values()
    
    mantel_r = mantelr(m1,m2)
    
    nrow = m1.shape[0]; ncol = m1.shape[1]
    
    # random permutation
    null_dist = np.zeros(args.nperm)
    
    tic = time.perf_counter()
    for i in range(args.nperm):
      permr = rng.permutation(nrow)
      permc = rng.permutation(ncol)
      null_dist[i] = mantelr(m1, m2[permr,:][:,permc])
      if i % 5000 == 0:
        toc = time.perf_counter() - tic
        print(f'{toc:.2f} seconds, {i}/{args.nperm}')
    
    # parse p value
    p = (null_dist > abs(mantel_r)).sum() + (null_dist < -abs(mantel_r)).sum()
    p /= args.nperm
      
    # print output
    print('Files correlated:')
    print(pf1, pf2, sep = '\n')
    print()
    print('Mantel R test statistic:')
    print(mantel_r)
    print('P-value:')
    print(p)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
      description = 'Extended Mantel test for phenotypic and genotypic correlation matrices '+
      '\nPlease supply two matrices of the same size and results will be printed to stdout')
    parser.add_argument('-i','--in', dest = '_in', help = 'input files', nargs = 2)
    parser.add_argument('--nperm', dest = 'nperm', help = '# permutations', type = int, 
                        default = 100000)
    args = parser.parse_args()
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    for i in args._in: proj.add_input(i,__file__)
    try: main(args)
    except: cmdhistory.errlog()