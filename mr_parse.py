#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2024-10-17

Parses Mendelian Randomisation for the same group of phenotypes.
Scans the entire directory for GWAS summary stats of the same data extension.
Creates summary tables for all MR results of different methods.
Significant targets identified from summary tables can be manually extracted for visual analysis
'''

def parse_mr_results(prefix):
    '''
    Data specifications:
        main_mr: 3 rows corresponding to MR IVW, WM and Egger
        presso: 3 rows corresponding to raw, outlier-corrected, and RSS-obs
        dirtest: one row, last two elements are 'correct direction' and Steiger p value
        pleio: one row, last three elements are 'egger intercept', SE and p value
    Outputs:
        causal output summary
        pleiotropy summary, by pleiotropy test and presso RSS-obs
    '''
    
    import pandas as pd
    
    main_mr = pd.read_table(f'{prefix}_results.txt')
    presso = pd.read_table(f'{prefix}_presso_results.txt')
    dirtest = pd.read_table(f'{prefix}_dirtest.txt')
    pleio = pd.read_table(f'{prefix}_pleiotropy.txt')
    
    # pleiotropy output: egger intercept and p value, rssobs and p value
    pleio = pleio[['outcome','exposure','egger_intercept','se','pval']]
    pleio.columns = ['outcome','exposure','egger_intercept','se','egger_p']
    pleio['rssobs'] = presso.iloc[-1,-2]
    pleio['presso_p'] = presso.iloc[-1,-1]
    
    # merge main and presso results to create the causal output summary table
    main_mr = main_mr.iloc[:,2:] # removes id.exposure and id.outcome
    presso.columns = ['method','b','se','t','pval']
    presso['outcome'] = main_mr['outcome']
    presso['exposure'] = main_mr['exposure']
    presso['nsnp'] = main_mr['nsnp']
    presso['method'] = ['MR-PRESSO raw','MR-PRESSO outlier-corrected','MR-PRESSO global']
    presso = presso.drop('t', axis = 'columns') #NB this step also drops the RSSobs estimate
    presso = presso[main_mr.columns]
    causal = pd.concat((main_mr, presso.iloc[:-1,:]), axis = 'index')
    causal['p'] = causal['pval']; causal = causal.drop('pval', axis = 1)
    
    causal['correct_dir'] = dirtest.iloc[0,-2]
    causal['dirtest_p'] = dirtest.iloc[0,-1]
    
    return causal, pleio

def stratified_fdr(df,label, pvalues):
    import pandas as pd
    import numpy as np
    from scipy.stats import false_discovery_control as fdr
    
    # df is a dataframe, label is a column that stratifies data, pvalues are column names that specify p-values
    labels = df[label].unique()
    groups_list = []
    grouped = df.groupby(label)
    
    for l in labels:
        group = grouped.get_group(l)
        for pcol in pvalues:
            p = group[pcol].to_numpy()
            p[p=='<0.001'] = '0' # MR-PRESSO outputs may output <0.001
            p = p.astype(np.float64)
            q = np.zeros(p.size)
            q[~pd.isna(group[pcol])] = fdr(p[~pd.isna(group[pcol])])
            q[pd.isna(group[pcol])] = np.nan
            qcol = pcol[:-1]+'pfdr'
            group[qcol] = q
        groups_list.append(group)
    out = pd.concat(groups_list)
    return out

def main(args):
    import os
    from fnmatch import fnmatch
    import pandas as pd
    
    for p2 in args.p2:
      # input processing by scanning directories
      prefix2 = []
      for f in os.listdir(f'{args.gwa}/{p2}'):
          if fnmatch(f, '*_X*') or fnmatch(f,'*_all_chrs*'): continue # exclude x-chr sum stats
          if fnmatch(f,f'*.{args.ext2}'): prefix2.append(f.replace(f'.{args.ext2}',''))
    
      for f2 in prefix2:
        for p1 in args.p1:
          # initialise output parsed tables
          results_fwd = []
          results_rev = []
          pleio_fwd = []
          pleio_rev = []    
        
          # input processing by scanning directories
          prefix1 = []
          for f in os.listdir(f'{args.gwa}/{p1}'):
              if fnmatch(f, '*_X*') or fnmatch(f,'*_all_chrs*'): continue # exclude x-chr sum stats
              if fnmatch(f, f'*.{args.ext1}'): prefix1.append(f.replace(f'.{args.ext1}',''))
          
          for f1 in prefix1:
            mr_prefix = f'{args._in}/{p2}/{f2}/{p1}_{f1}_{f2}'
            try: c, p = parse_mr_results(f'{mr_prefix}_mr_forward')
            except: print(f'{mr_prefix} no MR forward result')
            results_fwd.append(c); pleio_fwd.append(p)
            try: c, p = parse_mr_results(f'{mr_prefix}_mr_reverse')
            except: print(f'{mr_prefix} no MR reverse result')
            results_rev.append(c); pleio_rev.append(p)
        
          # concatenate and write tabular output
          results_fwd = pd.concat(results_fwd)
          results_fwd_corrected = stratified_fdr(results_fwd,'method',['p','dirtest_p'])
          results_fwd_corrected.sort_values(by = 'pfdr', inplace = True)
          results_fwd_corrected.to_csv(f'{args._in}/{p2}/{p1}_{f2}_mr_forward.txt', sep = '\t', index = False)
          
          results_rev = pd.concat(results_rev)
          results_rev_corrected = stratified_fdr(results_rev,'method',['p','dirtest_p'])
          results_rev_corrected.sort_values(by = 'pfdr', inplace = True)
          results_rev_corrected.to_csv(f'{args._in}/{p2}/{p1}_{f2}_mr_reverse.txt', sep = '\t', index = False)
          
          pleio_fwd = pd.concat(pleio_fwd)
          pleio_fwd_corrected = stratified_fdr(pleio_fwd,'outcome',['egger_p','presso_p'])
          pleio_fwd_corrected.to_csv(f'{args._in}/{p2}/{p1}_{f2}_mr_forward_pleiotropy.txt',
                                     sep = '\t', index = False)
          
          pleio_rev = pd.concat(pleio_rev)
          pleio_rev_corrected = stratified_fdr(pleio_rev,'exposure',['egger_p','presso_p'])
          pleio_rev_corrected.to_csv(f'{args._in}/{p2}/{p1}_{f2}_mr_reverse_pleiotropy.txt', 
                                     sep = '\t', index = False)
          
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'This script parses MR outputs for groups of phenotypes, please run mr_batch.py first')
    parser.add_argument('-g','--gwa', dest = 'gwa', 
                        help = 'input GWA directory, assumes both groups of pheno to be in the same dir',
                        default = '../gwa')
    parser.add_argument('-p1','--pheno1', dest = 'p1', 
                        help = 'Phenotypes group 1 (reserved for IDPs)', nargs = '*',
                        default=['deg_local','degi_local','degc_local',
                                 'clu_local','eff_local','mpl_local'])
    parser.add_argument('-e1','--ext1', dest = 'ext1', help = 'Extension for phenotype group 1',
                        default = 'fastGWA')
    parser.add_argument('-p2','--pheno2', dest = 'p2', 
                        help = 'Phenotypes group 2 (reserved for correlates)', nargs = '*', 
                        default = ['disorders_for_mr']) # require manual fiddling, so create new dir
    parser.add_argument('-e2','--ext2', dest = 'ext2', help = 'Extension for phenotype group 2',
                        default = 'txt')
    parser.add_argument('-i','--in', dest = '_in', help = 'Input directory, should be the output of mr_batch.py',
                        default = '../mr')
    # default output to the same directory
    parser.add_argument('-f','--force', dest = 'force', action = 'store_true',
                        default = False, help = 'Force overwrite')
    args = parser.parse_args()
    
    import os
    for arg in ['gwa','_in']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
    from _utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_input(f'{args._in}/{args.p2}/*/*',__file__)
    proj.add_output(f'{args._in}/{args.p2}/*.txt',__file__)
    
    try: main(args)
    except: cmdhistory.errlog()