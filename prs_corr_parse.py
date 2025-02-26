#!/usr/bin/env python

def main(args):
    import os
    import pandas as pd
    import scipy.stats as sts
    import numpy as np
    from _utils.path import normaliser
    norm = normaliser()
    
    summary = []
    for x in args.pheno:
        if not os.path.isfile(f'{args._in}/{x}_summary.txt'):
            print(f'{x} not found in prs correlation records, check')
            continue
        summary.append(pd.read_table(f'{args._in}/{x}_summary.txt'))
    summary = pd.concat(summary)
    norm.normalise(summary).to_csv(f'{args.out}/prscorr_'+'_'.join(args.pheno) + '.' + '_'.join(args.prs)+'.txt', sep = '\t', index = False)
    
    for prs in args.prs:
        prsdf = summary.loc[summary.pheno2 == prs,:]
        norm.normalise(prsdf).to_csv(f'{args.out}/{prs}_summary.txt', sep = '\t', index = False)
        prs_wide = prsdf.pivot_table(columns = 'group1', index = 'pheno1', values = 'beta')
        prs_wide.index.name = None; prs_wide.columns.name = None
        prs_wide.to_csv(f'{args.out}/{prs}_beta.txt', sep = '\t', index = False) # do not normalise to feed into brainplots
        
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
      description = 'this script prepares the PRS correlations for r-ggseg plotting')
    parser.add_argument('pheno', nargs = '*', help = 'local phenotypes',
      default=['deg_local','degi_local','degc_local','clu_local','eff_local','mpl_local'])
    parser.add_argument('-i','--in', dest = '_in', help = 'input directory',
      default='../prs/prs_corr/')
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory')
    parser.add_argument('--prs', dest = 'prs', help = 'PRS scores to gather', nargs = '*',
        default = ['an2019','asd2022','ptsd2024','mdd2023','scz2022','bip2021','adhd2022',
       'sud2023'])
    # always overwrites
    args = parser.parse_args()
    import os
    args._in = os.path.realpath(args._in)
    if type(args.out) == type(None): args.out = args._in
    args.pheno.sort()
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheno_.*.csv', __file__)
    proj.add_output(args.out+'/%pheno_.*.csv', __file__)
    try: main(args)
    except: cmdhistory.errlog()