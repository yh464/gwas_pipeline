#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2024-10-29

Summarises the fine-map results into a table of sig variants
'''

def main(args):
    import os
    import pandas as pd
    import numpy as np
    from fnmatch import fnmatch
    from _utils.path import normaliser
    norm = normaliser()
    
    dflist = []
    for x in args.pheno:
      os.chdir(args._in)
      os.chdir(x)
      
      flist = []
      for y in os.listdir():
        if fnmatch(y,'bet_asym*') or fnmatch(y,'*frac*') or fnmatch(y,'*asym_abs*'):
          continue
        if fnmatch(y,'*finemap.summary'):
          flist.append(y)
      
      for y in flist:
        prefix = y.split('.')[0].replace('_0','')
        df = pd.read_csv(y, sep = '\s+')
        df = df.loc[df.CREDIBLE_SET!=0,:]
        df.insert(loc = 0, column = 'PHENO', value = prefix)
        dflist.append(df)
    
      out_df = pd.concat(dflist, axis = 0).sort_values(by = ['CHR','BP'])
      norm.normalise(out_df).to_csv(f'{args.out}/{x}_sig_variants.txt', header = True, index = False, sep = '\t')
      np.savetxt(f'{args.out}/{x}_sig_variants_list.txt',out_df['SNP'].unique(),fmt = '%s')
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
      description = 'This script parses fine-mapping summary files')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*',
      default=['deg_local','degi_local','degc_local','clu_local','eff_local','mpl_local'])
    parser.add_argument('-i','--in', dest = '_in', help = 'Fine-mapping summary files directory', # in and out to same directory
      default = '../finemap/')
    parser.add_argument('-o', '--out', dest = 'out', help = 'output directory')
    # always overwrites
    args = parser.parse_args()
    # path normalisation
    
    import os
    for arg in ['_in',]:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    if type(args.out) == type(None): args.out = args._in
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheng/%pheno_%maf.finemap.summary',__file__)
    proj.add_output(args.out+'/sig_variants.csv',__file__)
    try: main(args)
    except: cmdhistory.errlog()