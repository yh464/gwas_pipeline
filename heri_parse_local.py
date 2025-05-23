#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2023-07-13
Version 2: 2024-11-14

Creates a heritability summary table for local phenotypes

Preceding workflow:
    heri_batch.py
Requires following inputs:
    LDSC H2 logs
'''

def main(args):
    import os
    from fnmatch import fnmatch
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    import scipy.stats as sts
    
    m2m = pd.read_csv('/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/params/hcp2yeo.csv')
    m2m = m2m[['label1','label2']]
    m2m.columns = ['roi','Yeo']
    
    all_summary = []
    for x in args.pheno:
        # if not fnmatch(x, '*local*'): continue
        os.chdir(args._in)
        glob_h2 = np.nan; glob_se = np.nan
        for y in os.listdir(args.glob):
            if fnmatch(y, x.replace('local','global')+'*.h2.log'):
                tmp = open(f'{args.glob}/{y}')
                for z in tmp.read().splitlines():
                    if fnmatch(z, 'Total Observed scale h2*'):
                        l = z.split(' ')
                        break
                try:
                    glob_h2 = float(l[-2])
                    glob_se = float(l[-1].replace('(','').replace(')',''))
                except:
                    glob_h2 = np.nan; glob_se = np.nan
        os.chdir(x)
        
        flist = []
        for y in os.listdir():
          if fnmatch(y, '*.h2.log'):
            flist.append(y)
        
        summary = []
        for y in flist:
          prefix = y.replace('.h2.log', '')
          name = prefix.replace('_0.01','')
          f = open(y)
          tmp = f.read().splitlines()
          
          for z in tmp:
            if fnmatch(z, 'Total Observed scale h2*'):
              l = z.split(' ')
              break
          try:
              h2 = float(l[-2])
              se = float(l[-1].replace('(','').replace(')',''))
          except:
              h2 = np.nan; se = np.nan
          summary.append(pd.DataFrame(dict(
              pheno = x, roi = [name], metric = 'heritability', method = 'LDSC', # placeholder to put h2 on 5th column
              h2 = h2, se = se, z = h2/se, p = 1-sts.chi2.cdf(h2**2/se**2, df = 1),
              glob_h2 = glob_h2, glob_se = glob_se)))
      
        s = pd.concat(summary)
        s['z'] = s.h2 / s.se
        all_summary.append(s)
        s = s.merge(m2m, on = 'roi', how = 'outer')
        ltmp = s.Yeo.values
        ltmp [ltmp == np.nan] = 'subcortical'
        s.to_csv(f'{args._in}/{x}/h2_summary.txt', sep = '\t', index = False)
        try:
            sns.displot(s, x = 'h2', col = 'Yeo', bins = 20)
            plt.savefig(f'{args._in}/{x}/h2_summary.png')
            plt.close()
            sns.histplot(s, x = 'z', bins = 20)
            plt.savefig(f'{args._in}/{x}/h2_summary_z.png')
            plt.close()
        except:
            print('Check naming conventions, no fig plotted')
    
    all_summary = pd.concat(all_summary)
    h2_all = all_summary.pivot_table(values = 'h2', index = 'roi', columns = 'pheno')
    h2_all.index.name = 'label'
    
    from _utils.path import normaliser
    norm = normaliser()
    norm.normalise(all_summary).to_csv(f'{args.out}/local_h2_summary.txt', sep = '\t', index = False)
    norm.normalise(h2_all).to_csv(f'{args.out}/local_h2_summary.wide.txt', sep = '\t', index = False)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'Creates a heritability summary table for local phenotypes')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*',
      default=['deg_local','degi_local','degc_local','clu_local','eff_local','mpl_local'])
    parser.add_argument('-g','--glob', dest = 'glob', help = 'global phenotype', default = 'global_graph')
    parser.add_argument('-i','--in', dest = '_in', help = 'LDSC h2 log directory, sub-directories are PHENO elements',
      default = '../gcorr/ldsc_sumstats/')
    parser.add_argument('-o','--out', dest = 'out', help = 'Summary table output directory',
      default = '../local_corr/')
    # always overwrites
    args = parser.parse_args()
    import os
    args._in = os.path.realpath(args._in)
    args.out = os.path.realpath(args.out)
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheng/%reg_%maf.h2.log', __file__)
    proj.add_output(args._in+'/local_h2_summary.txt', __file__)
    try: main(args)
    except: cmdhistory.errlog()