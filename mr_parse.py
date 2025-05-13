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
    import numpy as np
    
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
    try:
        cause = pd.read_table(f'{prefix}_cause_results.txt')
        causal['cause_p'] = cause.iloc[-1,-1]
    except: causal['cause_p'] = np.nan
    try:
        lcv = f'{prefix}_lcv_results.txt'.replace('_forward','').replace('_reverse','')
        lcv = pd.read_table(lcv)
        causal['lcv_p'] = lcv.loc['p','x']
    except: causal['lcv_p'] = np.nan
    
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
            q[~group[pcol].isna().to_numpy()] = fdr(p[~group[pcol].isna().to_numpy()])
            q[group[pcol].isna().to_numpy()] = np.nan
            qcol = pcol[:-1]+'q'
            q = pd.DataFrame(data = q, index = group.index, columns =[qcol])
            group = pd.concat([group, q],axis = 1)
        groups_list.append(group)
    out = pd.concat(groups_list)
    return out

def main(args):
    import warnings
    from fnmatch import fnmatch
    import pandas as pd
    from _utils.path import normaliser, find_gwas, pair_gwas
    
    norm = normaliser()
    exposures = find_gwas(args.p1, ext = args.ext1, se = True)
    outcomes = find_gwas(args.p2, ext = args.ext2, se = True)
    missing = []

    for g2, p2s in outcomes:
      for g1, p1s in exposures:
        all_fwd = []
        all_rev = []
        all_compare = []
        
        for p2 in p2s:
            # initialise output parsed tables
            results_fwd = []
            results_rev = []
            results_compare = []
            pleio_fwd = []
            pleio_rev = []    
            
            for p1 in p1s:
                mr_prefix = f'{args._in}/{g2}/{p2}/{g1}_{p1}_{p2}'
                try: 
                    cfwd, p = parse_mr_results(f'{mr_prefix}_mr_forward')
                    results_fwd.append(cfwd); pleio_fwd.append(p)
                    crev, p = parse_mr_results(f'{mr_prefix}_mr_reverse')
                    results_rev.append(crev); pleio_rev.append(p)
                    compare = pd.DataFrame(index = [0], 
                        columns = ['outcome','exposure','b_fwd','se_fwd','p_fwd','b_rev','se_rev','p_rev'])
                    compare.iloc[0,:5] = cfwd.loc[cfwd.method == 'MRlap corrected IVW',['outcome','exposure','b','se','p']].iloc[0,:]
                    compare.iloc[0,5:] = crev.loc[crev.method == 'MRlap corrected IVW',['b','se','p']].iloc[0,:]
                    results_compare.append(compare)
                except: missing.append(mr_prefix)
            
            # concatenate and write tabular output
            try:
                results_fwd = pd.concat(results_fwd)
                results_fwd_corrected = stratified_fdr(results_fwd,'method',['p','cause_p'])
                results_fwd_corrected.sort_values(by = 'q', inplace = True)
                all_fwd.append(results_fwd_corrected)
                results_fwd_corrected.to_csv(f'{args._in}/{g2}/{g1}_{p2}_mr_forward.txt', sep = '\t', index = False)
                pleio_fwd = pd.concat(pleio_fwd)
                pleio_fwd_corrected = stratified_fdr(pleio_fwd,'outcome',['egger_p','presso_p'])
                pleio_fwd_corrected.to_csv(f'{args._in}/{g2}/{g1}_{p2}_mr_forward_pleiotropy.txt',
                                            sep = '\t', index = False)
            except:
                warnings.warn(f'{p2} no MR forward results - check summary stats')
            
            try:
                results_rev = pd.concat(results_rev)
                results_rev_corrected = stratified_fdr(results_rev,'method',['p','cause_p'])
                results_rev_corrected.sort_values(by = 'q', inplace = True)
                all_rev.append(results_rev_corrected)
                results_rev_corrected.to_csv(f'{args._in}/{g2}/{g1}_{p2}_mr_reverse.txt', sep = '\t', index = False)
                pleio_rev = pd.concat(pleio_rev)
                pleio_rev_corrected = stratified_fdr(pleio_rev,'exposure',['egger_p','presso_p'])
                pleio_rev_corrected.to_csv(f'{args._in}/{g2}/{g1}_{p2}_mr_reverse_pleiotropy.txt', 
                                            sep = '\t', index = False)
            except:
                warnings.warn(f'{p2} no MR reverse results - check summary stats')
            
            try:
                results_compare = pd.concat(results_compare)
                results_compare.to_csv(f'{args._in}/{g2}/{g1}_{p2}_mr_compare.txt', sep = '\t', index = False)
                all_compare.append(results_compare)
            except: pass
        
        norm.normalise(pd.concat(all_fwd).sort_values(by = 'q')).to_csv(
            f'{args._in}/{g2}/all_{g1}_{g2}_mr_forward.txt', sep = '\t', index = False)
        norm.normalise(pd.concat(all_rev).sort_values(by = 'q')).to_csv(
            f'{args._in}/{g2}/all_{g1}_{g2}_mr_reverse.txt', sep = '\t', index = False)
        norm.normalise(pd.concat(all_compare)).to_csv(
            f'{args._in}/{g2}/all_{g1}_{g2}_mr_compare.txt', sep = '\t', index = False)
    print('Missing MR results:')
    for m in missing: print(m)

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
                        default = 'fastGWA')
    parser.add_argument('-i','--in', dest = '_in', help = 'Input directory, should be the output of mr_batch.py',
                        default = '../mr')
    # default output to the same directory
    # always overwrites
    args = parser.parse_args()
    args.p1.sort()
    args.p2.sort()
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