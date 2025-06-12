#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-01-23

Summarises summary-data Mendelian Randomisation outputs

Upstream workflow: 
    annot_smr_batch.py
'''

def read_smr(file):
    import pandas as pd
    df = pd.read_table(file, usecols = ['probeID','ProbeChr','Gene', 'topSNP','A1','A2','b_SMR','se_SMR','p_SMR','p_HEIDI','nsnp_HEIDI'])
    df.columns = ['probe','chr','gene','SNP','A1','A2','beta','se','p','p_heidi','nsnp_heidi']
    return df

def main(args):
    # parse input xqtl file
    from fnmatch import fnmatch
    import numpy as np
    import pandas as pd
    from scipy.stats import false_discovery_control as fdr
    from _utils.path import normaliser
    norm = normaliser()
    os.chdir(args._in)
    
    qtl_list = []
    for y in os.listdir(args.qtl):
        if fnmatch(y, '*.besd'): qtl_list.append(y)
        if os.path.isdir(f'{args.qtl}/{y}'):
            if any([fnmatch(z, '*.besd') for z in os.listdir(f'{args.qtl}/{y}')]):
                qtl_list.append(y)
    print('Following QTL have been found:')
    for qtl in qtl_list: print(qtl)
    print()    
    
    for x in args.pheno:
        # overall summary table
        all_phenos = []
        
        prefix_list = []
        for y in os.listdir(f'{args.gwa}/{x}'):
            if fnmatch(y,'*.fastGWA') and (not fnmatch(y,'*_X.fastGWA')) and (not fnmatch(y, '*_all_chrs*')):
                prefix_list.append(y.replace('.fastGWA',''))
        
        print('Following phenotypes have been found:')
        for y in prefix_list: print(y)
        
        for y in prefix_list:
            all_qtls = []
            for qtl in qtl_list:
                if os.path.isfile(f'{args._in}/{x}/{y}.{qtl}.smr'):
                    smr = read_smr(f'{args._in}/{x}/{y}.{qtl}.smr')
                else:
                    smr = []
                    for chrom in os.listdir(f'{args._in}/{x}/{y}.{qtl}'):
                        if not fnmatch(chrom, '*.smr'): continue
                        smr.append(read_smr(f'{args._in}/{x}/{y}.{qtl}/{chrom}'))
                    try: 
                        smr = pd.concat(smr)
                        smr.insert(loc = 0, column = 'qtl', value = qtl)
                        smr.insert(loc = 0, column = 'pheno', value = y)
                    except: 
                        smr = pd.DataFrame(columns = ['pheno','qtl','probe','chr','gene','SNP','A1','A2','beta',
                                                      'se','p','p_heidi','nsnp_heidi','q'], index = [])
                        all_qtls.append(smr)
                        print(f'WARNING: SMR results for {x}/{y}.{qtl} are missing')
                        continue
                smr['q'] = np.nan
                smr.loc[~smr.p.isna(),'q'] = fdr(smr.loc[~smr.p.isna(),'p'])
                all_qtls.append(smr)
            all_qtls = pd.concat(all_qtls).sort_values(by = ['q', 'p_heidi'])
            all_qtls.to_csv(f'{args._in}/{x}/{y}.txt', sep = '\t', index = False)
            all_phenos.append(all_qtls)
        all_phenos = pd.concat(all_phenos).sort_values(by = ['q','p_heidi'])
        all_phenos = norm.normalise(all_phenos)
        all_phenos.to_csv(f'{args._in}/{x}.txt', sep = '\t', index = False)
        all_phenos.loc[all_phenos.p < 0.05,:].to_csv(f'{args._in}/{x}_sig.txt', sep = '\t', index = False)
        
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description = 'This programme summarises summary data randomisation outputs')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
    parser.add_argument('-i','--in', dest = '_in', help = 'Directory containing all SMR summary statistics',
        default = '../annot/smr')
    parser.add_argument('--gwa', dest = 'gwa', 
        help = 'Directory containing all GWA summary statistics, just for scanning phenotypes',
        default = '../gwa/')
    parser.add_argument('-q','--qtl', dest = 'qtl', help = 'Directory containing all xQTL files',
        default = '../params/xqtl')
    # always overwrites
    args = parser.parse_args()
    import os
    for arg in ['_in','qtl','gwa']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    
    from _utils import cmdhistory, logger
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()