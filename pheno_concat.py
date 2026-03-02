#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1.0: 2024-11-26

A general utility to concatenate phenotype files of all subjects
due to different phenotype file formats, a pheno_concat_legacy file is preserved

Preceding workflow:
    pheno_batch.py (and other phenotyping scripts)
Requires following inputs: 
    phenotype file for each subject in following format:
        tab-separated
        index = phenotypes (to be included as columns in the output files)
        columns (header) = phenotype groups (to be file names of the output files)
Output:
    file name: phenotype group; index = subject ID; columns = phenotypes
    downstream: 

Changelog:
    changed input format so that index = phenotype name, columns = phenotype group name
'''

import pandas as pd
import numpy as np
import os
from _utils.logger import logger
log = logger()

def read_file(file):
    subj = os.path.basename(file).replace('.txt','')
    df = pd.read_table(file, index_col = 0)
    df['pheno'] = df.index
    df = df.melt(id_vars = 'pheno', var_name = 'pheng')
    df.insert(0, column = 'EID', value = subj)
    return df

def main(args):
    from tqdm import tqdm
    from multiprocessing import Pool
    from fnmatch import fnmatch
    os.chdir(args._in)
    pool = Pool(64)
    for x in args.pheno:
        files = [f'{x}/{y}' for y in os.listdir(x) if fnmatch(y, '*.txt')]
        dflist = list(tqdm(pool.imap(read_file, files, chunksize = 64), total = len(files), desc = f'Processing {x}'))
        df = pd.concat(dflist)
        pheno_groups = df['pheng'].unique()
        for pg in pheno_groups:
            tmp = df.loc[df.pheng == pg, :]
            tmp = tmp.pivot_table(columns = 'pheno', index = 'EID', values = 'value')
            tmp.to_csv(f'{pg}.unstd.txt', sep = '\t', index = False)
            for c in tmp.columns:
                if (tmp[c] == np.inf).any():
                    log.warn(f'Phenotype {pg}/{c} has infinite values, setting to NaN')
                    tmp[c] = tmp[c].replace(np.inf, np.nan)
                if tmp[c].std() <= 0: log.warn(f'Phenotype {pg}/{c} has zero variance and will not be scaled'); continue
            tmp = tmp.copy()
            tmp /= tmp.std(axis = 0)
            tmp.insert(0, column = 'FID', value = tmp.index)
            tmp.insert(1, column = 'IID', value = tmp.index)
            tmp.to_csv(f'{pg}.txt', index = False, sep = '\t')
    return

if __name__ == '__main__':
    # input argument processing
    import argparse as ap
    parser = ap.ArgumentParser(description='This programme concatenates all subjects for their functional connectome phenotypes')
    parser.add_argument('pheno', nargs = '*', help = 'Phenotypes to concatenate')
    parser.add_argument('-i','--in', dest = '_in', 
                        default = '../pheno/ukb/',
                        help = 'Data directory')
    parser.add_argument('-f','--force', dest = 'force', help = 'Force output',
                        default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    args._in = os.path.realpath(args._in)
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in, __file__)
    proj.add_output(args._in, __file__)
    try: main(args)
    except: cmdhistory.errlog()