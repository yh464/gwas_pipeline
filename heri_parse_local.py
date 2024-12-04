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
    
    m2m = pd.read_csv('/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/params/map2map.csv')
    m2m = m2m[['label1','label2']]
    m2m.columns = ['roi','Yeo']
    
    all_h2 = []
    all_z = []
    for x in args.pheno:
      if not fnmatch(x, '*local*'): continue
      os.chdir(args._in)
      os.chdir(x)
      
      flist = []
      for y in os.listdir():
        if fnmatch(y, '*.h2.log'):
          flist.append(y)
      
      namelist = []
      h2list = []
      selist = []
      for y in flist:
        prefix = y.replace('.h2.log', '')
        namelist.append(prefix.replace('_0.01',''))
        f = open(y)
        tmp = f.readlines()
        
        for z in tmp:
          if fnmatch(z, 'Total Observed scale h2*'):
            l = z.split(' ')
            break
        
        h2list.append(float(l[-2]))
        se = l[-1].replace('(','').replace(')','')
        selist.append(float(se))
    
      s = pd.DataFrame(dict(roi = namelist, h2 = h2list, se = selist))
      s['z'] = s.h2 / s.se
      all_h2.append(s.h2.rename(x))
      all_z.append(s.z.rename(x))
      s = s.merge(m2m, on = 'roi', how = 'outer')
      ltmp = s.Yeo.values
      ltmp [ltmp == np.nan] = 'subcortical'
      s.to_csv(f'{args._in}/{x}/h2_summary.txt', sep = '\t', index = False)
      sns.displot(s, x = 'h2', col = 'Yeo', bins = 20)
      plt.savefig(f'{args._in}/{x}/h2_summary.png')
      plt.close()
      sns.histplot(s, x = 'z', bins = 20)
      plt.savefig(f'{args._in}/{x}/h2_summary_z.png')
      plt.close()
    
    tmp = namelist
    glasser = []
    for i in tmp:
      i = i.replace('_ROI','')
      i = 'lh_'+ i if i[0]=='L' else 'rh_'+i
      glasser.append(i)
    glasser = pd.Series(glasser, name = 'label')
    h2_all = pd.concat([glasser]+all_h2, axis = 1).sort_values(by = 'label')
    h2_all.to_csv(f'{args.out}/local_h2_summary.csv', sep = ',', index = False)
    z_all = pd.concat([glasser]+all_z, axis = 1).sort_values(by = 'label')
    z_all.to_csv(f'{args.out}/local_z_summary.csv', sep = ',', index = False)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'Creates a heritability summary table for local phenotypes')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*',
      default=['deg_local','degi_local','degc_local','clu_local','eff_local','mpl_local'])
    parser.add_argument('-i','--in', dest = '_in', help = 'LDSC h2 log directory, sub-directories are PHENO elements',
      default = '../gene_corr/ldsc_sumstats/')
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