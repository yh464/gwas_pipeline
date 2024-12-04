#!/usr/bin/env python

def main(args):
    import os
    import pandas as pd
    
    os.chdir(args._in)
    
    for prs in args.prs:
      exec(f'beta_{prs} = []; z_{prs} = []')
    
    for x in args.pheno:
      beta_df = pd.read_csv(f'{x}_beta.csv', sep = '\s+', index_col = 0)
      for prs in args.prs:
        exec(f'beta_{prs}.append(beta_df[prs])')
      
      z_df = pd.read_csv(f'{x}_z.csv', sep = '\s+', index_col = 0)
      for prs in args.prs:
        exec(f'z_{prs}.append(z_df[prs])')
    
    import scipy.stats as sts
    import numpy as np
    
    for prs in args.prs:
      tmp = eval(f"pd.concat(beta_{prs}, axis = 1)")
      tmp.index.name = 'label'
      tmp.columns = args.pheno
      tmp.to_csv(f'{args.out}/{prs}_beta.csv')
      tmp = eval(f"pd.concat(z_{prs}, axis = 1)")
      tmp.index.name = 'label'
      tmp.columns = args.pheno
      tmp.to_csv(f'{args.out}/{prs}_z.csv')
      p = 1 - sts.chi2.cdf(df = 1, x = tmp**2)
      p[np.isnan(p)] = 1
      pfdr = sts.false_discovery_control(p)
      pfdr[pfdr == 1] = np.nan
      pfdr = pd.DataFrame(pfdr, columns = tmp.columns, index = tmp.index)
      pfdr.to_csv(f'{args.out}/{prs}_pfdr.csv')
      
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
      description = 'this script prepares the PRS correlations for r-ggseg plotting')
    parser.add_argument('pheno', nargs = '*', help = 'local phenotypes',
      default=['deg_local','degi_local','degc_local','clu_local','eff_local','mpl_local'])
    parser.add_argument('-i','--in', dest = '_in', help = 'input directory',
      default='../prs/prs_corr/')
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory')
    parser.add_argument('--prs', dest = 'prs', help = 'PRS scores to gather, separate with comma',
      default = 'an2019,asd2019,ptsd2024,mdd2023,scz2022,bip2021,adhd2022,sud2023')
    # always overwrites
    args = parser.parse_args()
    import os
    args._in = os.path.realpath(args._in)
    if type(args.out) == type(None): args.out = args._in
    args.prs = args.prs.split(',')
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheno_.*.csv', __file__)
    proj.add_output(args.out+'/%pheno_.*.csv', __file__)
    try: main(args)
    except: cmdhistory.errlog()