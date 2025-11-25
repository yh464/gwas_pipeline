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
def main(args):
    import pandas as pd
    import os
    from fnmatch import fnmatch
    os.chdir(args._in)
    for x in args.pheno:
        dflist = []
        for y in os.listdir(x):
            if not fnmatch(y, '*.txt'): continue
            subj = y.replace('.txt','')
            df = pd.read_table(f'{x}/{y}', index_col = 0)
            df['pheno'] = df.index
            df = df.melt(id_vars = 'pheno', var_name = 'pheng')
            df.insert(0, column = 'EID', value = subj)
            dflist.append(df)
        df = pd.concat(dflist)
        pheno_groups = df['pheng'].unique()
        for pg in pheno_groups:
            tmp = df.loc[df.pheng == pg, :]
            tmp = tmp.pivot_table(columns = 'pheno', index = 'EID', values = 'value')
            for c in tmp.columns:
                tmp[c] /= tmp[c].std()
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
    
    from ._utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in, __file__)
    proj.add_output(args._in, __file__)
    try: main(args)
    except: cmdhistory.errlog()