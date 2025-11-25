#!/usr/bin/env python

def decomp(corrmat):
    import numpy as np
    corrmat = np.abs(corrmat)
    l,_ = np.linalg.eig(corrmat)
    n = l.size
    n_old = n * (1 - (n-1) * np.var(l) / n**2)
    p_old = 1-0.95**(1/n_old)
    l[l<0] = 0
    n_new = n * (1 - (n-1) * np.var(l) / n**2)
    p_new = 1-0.95**(1/n_new)
    return [n_old, p_old, n_new, p_new]

def main(args):
    import os
    import pandas as pd
    from fnmatch import fnmatch
    import numpy as np
    from time import perf_counter as t
    
    tic = t()
    # scan directory
    flist = []
    for p in args.pheno:
        # exact match
        if os.path.isfile(f'{args._in}/{p}.txt'):
            flist.append(f'{args._in}/{p}.txt')
            continue
        elif os.path.isfile(f'{args._in}/{p}'):
            flist.append(f'{args._in}/{p}')
            continue
        
        # else scan directory
        for x in os.listdir(args._in):
            if fnmatch(x, f'*{p}*.txt'):
                flist.append(f'{args._in}/{x}')
    print('Following phenotype files have been found:')
    for f in flist: print(f)
    
    pflist = []; dflist = []
    for x in flist:
        pflist.append('.'.join(os.path.basename(x).split('.')[:-1]))
        df = pd.read_table(x, sep = '\s+', index_col=['FID','IID'])
        for col in df.columns:
            if sum(df[col].isna()) > 0.5 * df.shape[0]:
                df.drop(col, axis = 1, inplace = True)
        dflist.append(df)
    
    dflist = pd.Series(dflist)
    n_pheno = len(dflist)
    
    if args._all:
        comb = []
        from itertools import combinations
        for m in range(1,n_pheno + 1):
            for c in combinations(range(n_pheno),m):
                tmp = np.zeros(n_pheno).astype('?')
                tmp[list(c)] = True
                comb.append(tmp.copy())
    else:
        comb = np.vstack((np.identity(n_pheno), np.ones((1,n_pheno)))).astype('?')
    
    comb_df = pd.DataFrame(data = comb, columns = pflist)
    out_df = pd.DataFrame(data = float(0), index = comb_df.index, 
                     columns = ['n_old','p_old','n_new','p_new'])
    
    for i in comb_df.index:
        c = pd.concat(dflist[comb_df.loc[i,:].tolist()].tolist(), axis = 1).corr()
        c = np.abs(c)
        l,_ = np.linalg.eig(c)
        n = l.size
        neff_abs = n * (1 - (n-1) * np.var(l) / n**2)
        out_df.loc[i,'n_old'] = neff_abs
        out_df.loc[i,'p_old'] = 1-0.95**(1/neff_abs)
        l[l<0] = 0
        neff_pos = n * (1 - (n-1) * np.var(l) / n**2)
        out_df.loc[i,'n_new'] = neff_pos
        out_df.loc[i,'p_new'] = 1-0.95**(1/neff_pos)
        toc = t() - tic
        print(f'{i+1}/{comb_df.shape[0]}, time = {toc:.3f}')
        
    out_df = pd.concat([comb_df, out_df], axis = 1)
    out_df.to_csv(f'{args._in}/neff_'+'_'.join(args.pheno)+'.txt', sep = '\t', index   = False, header = True)
    
    # with open(args.out,'w') as f:
    #   print(f'Effective # variables by old method Nyholt DR (2004): {neff_abs}', file = f)
    #   print(f'Bonferroni corrected threshold: {1-0.95**(1/neff_abs)}', file = f)
    #   print(file = f)
    #   print(f'Effective # variables by new method Nyholt DR (2004): {neff_pos}', file = f)
    #   print(f'Bonferroni corrected threshold: {1-0.95**(1/neff_pos)}', file = f)
      
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
      description = 'This programme calculates the effective number of indep. variables.')
    parser.add_argument('pheno', nargs = '*', help = 'Phenotypes to scan for')
    parser.add_argument('-i','--in', dest = '_in', help = 'List containing all phenotypes',
      default = '../pheno/ukb/')
    parser.add_argument('-a', dest = '_all', help = 'Iterate over all combinations of phenotypes',
      default = False, action = 'store_true')
    # always overwrites
    args = parser.parse_args()
    import os
    args._in = os.path.realpath(args._in)
    
    from ._utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in, __file__)
    proj.add_output(f'{args._in}/neff.txt', __file__)
    try: main(args)
    except: cmdhistory.errlog()