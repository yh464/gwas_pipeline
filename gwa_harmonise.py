#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-01-07
Version 2: 2025-04-23

Harmonises GWAS to fill in missing info from UKB genetics data

Requires following inputs: 
    GWAS summary statistics (scans directory)
    require the following columns: SNP, A1, A2, BETA or OR (manually edit if necessary)
'''

def harmonise(file, ref):
    import pandas as pd
    from fnmatch import fnmatch
    from time import perf_counter
    import warnings
    if type(ref) == str: ref = pd.read_table(ref, index_col='SNP')
    
    tic = perf_counter()
    if not fnmatch(file, '*.fastGWA') and not fnmatch(file,'*.txt'): return
    # read df
    try:
        df = pd.read_table(file, sep = '\\s+', index_col = 'SNP')
        toc = perf_counter() - tic
        print(f'Read {file}, time = {toc:.2f} seconds')
        to_drop = ['CHR','POS','BP']
        for col in ['CHR','POS','BP']: # erase CHR and BP information in case of different GRCh builds
            if not col in df.columns: to_drop.remove(col)
        df = df.drop(to_drop, axis = 1)
        for col in ['A1','A2']:
            df[col] = df[col].str.upper()
        df = df.rename(columns = {'A1':'A1_in','A2':'A2_in'})
        toc = perf_counter() - tic; print(f'Pre-processed {file}, time = {toc:.2f} seconds')

        # pre-process reference snp info
        ref = ref[['CHR','POS','A1','A2','AF1']] if not 'AF1' in ref.columns else ref[['CHR','POS','A1','A2']]
        
        df = pd.concat([ref, df], axis = 1, join = 'inner')
        matched = (df.A1 == df.A1_in) & (df.A2 == df.A2_in)
        flipped = (df.A1 == df.A2_in) & (df.A2 == df.A1_in)
        mismatch = (~matched) & (~flipped)
        toc = perf_counter() - tic; print(f'Identified mismatch for {file}, time = {toc:.2f} seconds')
        
        if 'OR' in df.columns: df.loc[flipped,'OR'] = df.loc[flipped,'OR'] ** -1
        if 'BETA' in df.columns: df.loc[flipped,'BETA'] *= -1
        if 'AF1' in df.columns: df.loc[flipped,'AF1'] = 1 - df.loc[flipped,'AF1']
        df = df.loc[~mismatch,:].drop(['A1_in','A2_in'], axis = 1).sort_values(by = ['CHR','POS'])
        toc = perf_counter() - tic
        print(f'Harmonised {file}, time = {toc:.2f} seconds')

        df.to_csv(file, sep = '\t', index = False)
        toc = perf_counter() - tic
        print(f'Saved {file}, time = {toc:.2f} seconds')
    except:
        warnings.warn(f'Failed to harmonise {file}, consider manual harmonisation')
    return

def main(args):
    import os
    import numpy as np
    import pandas as pd
    from fnmatch import fnmatch
    from multiprocessing import Pool
    from functools import partial
    
    # reference SNP info: CHR, SNP, POS, A1, A2, AF1
    fields = ['CHR','SNP','POS','A1','A2','AF1']
    ref = pd.read_table(args.ref, index_col='SNP')
    print('Read reference')
    flist = []
    for p in args.pheno:
        for file in os.listdir(f'{args._in}/{p}'):
            if fnmatch(f'{args._in}/{p}/{file}', '*.fastGWA') or fnmatch(file,'*.txt'):
                l = open(f'{args._in}/{p}/{file}').readline()
                if not 'SNP' in l or (not 'BETA' in l and not 'OR' in l): continue
                # check progress
                if all([x in l for x in fields]) and not args.force: continue
                flist.append(f'{args._in}/{p}/{file}')
    if len(flist) == 0: return
    pool = Pool(min((len(flist)),4))
    pool.map(partial(harmonise, ref=ref), flist, chunksize = int(np.ceil(len(flist)/4)))

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