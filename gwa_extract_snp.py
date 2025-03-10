#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-01-20

Utility script to extract input SNPs from GWAS files

Requires following inputs: 
    GWAS summary statistics (single file)
'''

def main(args):
    import os
    from fnmatch import fnmatch
    import pandas as pd
    
    all_files = []
    
    for x in args.pheno:
        os.chdir(args._in)
        os.chdir(x)
      
        # filter out required GWA files
        for y in os.listdir():
            if fnmatch(y, '*_X.fastGWA'): continue
            if fnmatch(y, '*_all_chrs.fastGWA'): continue
            if not fnmatch(y, '*.fastGWA') or fnmatch(y, '*.txt') or fnmatch(y, '*.tsv'): continue
            df = pd.read_table(y, sep = '\s+', usecols = ['SNP','BETA','SE','P'], index_col = 'SNP')
            df = df.loc[[args.snp],:].reset_index(names='SNP')
            prefix = '.'.join(y.split('.')[:-1]) # remove extension
            df['Phenotype'] = prefix
            df['Group'] = x
            df['Z'] = df.BETA/df.SE
            all_files.append(df[['Group','Phenotype','SNP','BETA','Z','P']])
        
    all_files = pd.concat(all_files)
    if args.out != None: all_files.to_csv(args.out, sep = '\t', index = False)
    all_files.to_clipboard(index = False)
    all_files_wide = all_files.pivot_table(values = 'P', index = 'SNP', columns = 'Phenotype')
    all_files_wide.insert(loc = 0, column = 'min_pval', value = all_files.min(axis = 1))
    print(all_files_wide)
    
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'This script generates PRS by continuous shrinkage from external sumstats')
    parser.add_argument('snp',help = 'list of SNPs to extract', nargs = '*')
    parser.add_argument('-p', dest = 'pheno', help = 'list of phenotype groups to extract SNPs',
      nargs = '*', default = ['global','global_asym'])
    parser.add_argument('-i','--in', dest = '_in', help = 'Directory containing all GWA summary statistics',
      default = '../gwa/')
    parser.add_argument('-o','--out', dest = 'out', help = 'Output list file')
    args = parser.parse_args()
    import os
    args._in = os.path.realpath(args._in)
    
    from _utils import cmdhistory, logger
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()