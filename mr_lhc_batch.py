#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2024-11-05

Batch submits jobs for Mendelian Randomisation for all GWAS files in a directory
(usually the same group of phenotypes). Scans the entire directory for GWAS summary
stats of the same data extension.
This script is specific to MR by LCV because it uses a different set of inputs
'''

def main(args):
    import os
    from fnmatch import fnmatch
    import numpy as np
    import pandas as pd
    
    # array submitter
    from _utils import array_submitter
    submitter = array_submitter.array_submitter(
        name = 'mr_lhc_'+'_'.join(args.p2), env = 'gentoolsr',
        timeout = 240, n_cpu = 2,
        debug = False
        )
    
    # output directory
    if not os.path.isdir(args.out): os.mkdir(args.out)
    
    force = '-f' if args.force else ''
    
    for p2 in args.p2:
      # make output directory wrt p2
      if not os.path.isdir(f'{args.out}/{p2}'): os.mkdir(f'{args.out}/{p2}')
        
      # input processing by scanning directories
      prefix2 = []
      for f in os.listdir(f'{args._in}/{p2}'):
          if fnmatch(f, '*_X*') or fnmatch(f,'*_all_chrs*'): continue # exclude x-chr sum stats
          if fnmatch(f,f'*.{args.ext2}'): prefix2.append(f.replace(f'.{args.ext2}',''))
      
      # parse n2 table
      if os.path.isfile(args.n2):
          n2_tbl = pd.read_table(args.n2)
          if n2_tbl.shape[1] == 2:
              n2_tbl.columns = ['pheno','n']
              n2_tbl['nca'] = np.nan
              n2_tbl['nco'] = np.nan
          elif n2_tbl.shape[1] == 4:
              n2_tbl.columns = ['pheno','n','nca','nco']
          else: raise ValueError('N2 should specify either both NCa and NCo or neither')
      else:
          if type(args.nca) != type(None):
              nca = args.nca; nco = args.n2 - args.nca
          else: nca = np.nan; nco = np.nan
          n2_tbl = pd.DataFrame(dict(pheno = prefix2, n = args.n2, nca = nca, nco = nco))
      n2_tbl.index = n2_tbl.pheno
      
      for p1 in args.p1:
        # input processing by scanning directories
        prefix1 = []
        for f in os.listdir(f'{args._in}/{p1}'):
            if fnmatch(f, '*_X*') or fnmatch(f,'*_all_chrs*'): continue # exclude x-chr sum stats
            if fnmatch(f, f'*.{args.ext1}'): prefix1.append(f.replace(f'.{args.ext1}',''))
        
        # parse n1 table
        if os.path.isfile(args.n1):
            n1_tbl = pd.read_table(args.n1)
            n1_tbl.columns = ['pheno','n']
        else:
            n1_tbl = pd.DataFrame(dict(pheno = prefix1, n = args.n1))
        n1_tbl.index = n1_tbl.pheno
        
        for f1 in prefix1:
          for f2 in prefix2:
            if not os.path.isdir(f'{args.out}/{p2}/{f2}'): os.mkdir(f'{args.out}/{p2}/{f2}')  
            
            # sample size
            n2 = n2_tbl.loc[f2,'n']
            n1 = n1_tbl.loc[f1, 'n']
            
            # check progress
            out_prefix = f'{args.out}/{p2}/{f2}/{p1}_{f1}_{f2}'
            
            # input file name specification
            gwa1 = f'{args._in}/{p1}/{f1}.{args.ext1}'
            gwa2 = f'{args._in}/{p2}/{f2}.{args.ext2}'
            
            if not os.path.isfile(f'{out_prefix}_mr_lhc_results.txt') or args.force:
                submitter.add(
                    f'Rscript mr_lhc.r --g1 {gwa1} --n1 {n1} --g2 {gwa2} --n2 {n2} '+
                    f'--ldsc {args.ldsc} --rho {args.rho} --hm3 {args.hm3} --refld {args.refld} '
                    f'-o {args.out}/{p2}/{f2} {force}')

    submitter.submit()
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'This script batch runs MR for groups of phenotypes')
    
    # input GWA
    g1 = parser.add_argument_group('GWAS parameters')
    g1.add_argument('-i','--in', dest = '_in', 
                    help = 'input GWA directory, assumes both groups of pheno to be in the same dir',
                    default = '../gwa')
    g1.add_argument('-p1','--pheno1', dest = 'p1', 
                    help = 'Phenotypes group 1 (reserved for IDPs)', nargs = '*',
                    default=['deg_local','degi_local','degc_local',
                             'clu_local','eff_local','mpl_local'])
    g1.add_argument('-e1','--ext1', dest = 'ext1', help = 'Extension for phenotype group 1',
                    default = 'fastGWA') # intended to scan the input directory
    g1.add_argument('-n1', help = 'Sample size for phenotype group 1, number or tabular',
                    default = 54030)
    g1.add_argument('-p2','--pheno2', dest = 'p2', 
                    help = 'Phenotypes group 2 (reserved for correlates)', nargs = '*', 
                    default = ['disorders_for_mr']) # require manual fiddling, so create new dir
    g1.add_argument('-e2','--ext2', dest = 'ext2', help = 'Extension for phenotype group 2',
                    default = 'txt')
    g1.add_argument('-n2', help = 'Sample size for phenotype group 2, number or tabular',
                    default = '../params/disorder_sample_size.txt')
    
    # LHC specific parameters
    g2 = parser.add_argument_group('LDSC parameters')
    g2.add_argument('--ldsc', help = 'LD score file (L2), better independent from study cohorts',
                    default = '../params/ldsc_for_gsem/uk10k.l2.ldscore')
    g2.add_argument('--rho',help = 'local LD score file (rho), better independent from study cohorts',
                    default = '../params/ldsc_for_gsem/uk10k.rho.ldscore')
    g2.add_argument('--hm3', help = 'HapMap3 SNP list',
                    default = '../params/ldsc_for_gsem/w_hm3.snplist')
    g2.add_argument('--refld', help = 'reference LD files for ancestry group, directory of 22 gz files',
                    default = '/rds/user/yh464/hpc-work/ldsc/baseline/') # intentionally absolute
    
    # output and force
    parser.add_argument('-o','--out', dest = 'out', help = 'Output directory',
                        default = '../mr')
    parser.add_argument('-f','--force', dest = 'force', action = 'store_true',
                        default = False, help = 'Force overwrite')
    args = parser.parse_args()
    
    import os
    for arg in ['_in','out','ldsc','rho','hm3','refld']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    try: args.n1 = int(args.n1)
    except: args.n1 = os.path.realpath(args.n1)
    try: args.n2 = int(args.n2)
    except: args.n2 = os.path.realpath(args.n2)
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(f'{args._in}/{args.p1}/*.{args.ext1}', __file__)
    proj.add_input(f'{args._in}/{args.p2}/*.{args.ext2}', __file__)
    proj.add_output(f'{args.out}/{args.p2}/*',__file__)
    
    try: main(args)
    except: cmdhistory.errlog()