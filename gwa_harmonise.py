#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-01-07

Harmonises GWAS to fill in missing info from UKB genetics data

Requires following inputs: 
    GWAS summary statistics (scans directory)
    require the following columns: SNP, A1, A2, BETA or OR (manually edit if necessary)
'''

def harmonise(file, ref):
    import pandas as pd
    from fnmatch import fnmatch
    if type(ref) == str: ref = pd.read_table(ref)
    
    if not fnmatch(file, '*.fastGWA') or fnmatch(file,'*.txt'): return
    # read df
    df = pd.read_table(file, sep = '\s+')
    for col in ['CHR','POS','BP']: # erase CHR and BP information in case of different GRCh builds
        if col in df.columns: df.drop(col, axis = 1, inplace = True)
    for col in ['A1','A2']:
        df[col] = df[col].str.upper()
        
    # construct df with flipped alleles
    df_rev = df.copy()
    df_rev.loc[:,'A2'] = df.loc[:,'A1'].copy() # intentional
    df_rev.loc[:,'A1'] = df.loc[:,'A2'].copy()
    if 'OR' in df.columns: df_rev['OR'] = df_rev['OR'] ** -1
    if 'BETA' in df.columns: df_rev['BETA'] *= -1
    
    if 'AF1' in df.columns:
        tmpref = ref.drop('AF1', axis = 1)
        df_rev.loc[:,'AF1'] *= -1
        df_rev.loc[:,'AF1'] += 1
    else: tmpref = ref
    
    merge = pd.merge(tmpref, df, on = ['SNP','A1','A2'])
    merge_rev = pd.merge(tmpref, df_rev, on = ['SNP','A1','A2'])
    
    df = pd.concat([merge, merge_rev], axis = 0).sort_values(by = ['CHR','POS'])
    df.to_csv(file, sep = '\t', index = False)
    

def main(args):
    import os
    import numpy as np
    import pandas as pd
    from fnmatch import fnmatch
    from multiprocessing import Pool
    from functools import partial
    
    # reference SNP info: CHR, SNP, POS, A1, A2, AF1
    fields = ['CHR','SNP','POS','A1','A2','AF1']
    ref = pd.read_table(args.ref)
    print('Read reference')
    for p in args.pheno:
        print(p)
        flist = []
        for file in os.listdir(f'{args._in}/{p}'):
            if fnmatch(f'{args._in}/{p}/{file}', '*.fastGWA') or fnmatch(file,'*.txt'):
                print(' '*4+file)
                l = open(f'{args._in}/{p}/{file}').readline()
                if not 'BETA' in l and not 'OR' in l: continue
                # check progress
                if all([x in l for x in fields]) and not args.force: continue
                flist.append(f'{args._in}/{p}/{file}')
                
        if len(flist) == 0: continue
        pool = Pool(min((len(flist)),10))
        pool.map(partial(harmonise, ref=ref), flist, chunksize = int(np.ceil(len(flist)/10)))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='This script harmonises GWAS summary stats to align with UKB genetics data')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
    parser.add_argument('-i','--in', dest = '_in', help = 'Input directory',
        default = '../gwa/')
    parser.add_argument('-r','--ref', dest = 'ref', help = 'reference UKB genetics data',
        default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/ukb_snp_info.txt')
    parser.add_argument('-f','--force', dest = 'force', action = 'store_true',
        default = False, help = 'force overwrite')
    args = parser.parse_args()
    import os
    args._in = os.path.realpath(args._in)
    
    from _utils import cmdhistory, logger
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()