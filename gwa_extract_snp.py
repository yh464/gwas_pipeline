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
    import pandas as pd
    from _utils.path import normaliser
    
    norm = normaliser()
    all_files = []
    tmpdir = os.path.realpath('../temp')
    
    for x in args.pheno:
        os.chdir(args._in)
        os.chdir(x)
        print(x)
        # search for required SNP using grep
        hdr = open(os.listdir()[0]).readline().replace('\n','').split()
        hdr[0] = 'Phenotype'
        for snp in args.snp:
            cache = f'{tmpdir}/sigsnp_{snp}_{x}.txt'
            if os.path.isfile(cache):
                df = pd.read_table(cache)
                all_files.append(df)
                continue
            os.system(f'grep {snp} *.fastGWA | xclip -selection clipboard') # quick search for selected SNP
            df = pd.read_clipboard(header = None)
            df.columns = hdr
            print(df.head())
            df = df.loc[df.SNP == snp,['Phenotype','SNP','BETA','SE','P']]
            df.columns = ['Phenotype','SNP','beta','se','p']
            phen = df.Phenotype
            phen = [y.split('.fastGWA')[0] for y in phen]
            df['Phenotype'] = phen
            df = df.loc[
                (~df.Phenotype.str.endswith('all_chrs')) & 
                (~df.Phenotype.str.endswith('_X')),:]
            df.insert(loc = 0, column = 'Group', value = x)
            df.insert(loc = 2, column = 'q', value = df.p * 1e+6)
            df.insert(loc = 4, column = 'Z', value = df.beta/df.se)
            df.to_csv(cache, sep = '\t', index = False)
            all_files.append(df)
            
    all_files = pd.concat(all_files)
    if args.out != None: norm.normalise(all_files).to_csv(args.out, sep = '\t', index = False)
    norm.normalise(all_files).to_clipboard(index = False)
    all_files_wide = all_files.pivot_table(values = 'p', index = 'snp', columns = 'Phenotype')
    all_files_wide.insert(loc = 0, column = 'min_pval', value = all_files.min(axis = 1))
    print(norm.normalise(all_files_wide))
    
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'This script extract selected SNPs from GWAS sumstats')
    parser.add_argument('snp',help = 'list of SNPs to extract', nargs = '*')
    parser.add_argument('-p', dest = 'pheno', help = 'list of phenotype groups to extract SNPs',
      nargs = '*', default = ['global','global_asym'])
    parser.add_argument('-i','--in', dest = '_in', help = 'Directory containing all GWA summary statistics',
      default = '../gwa/')
    parser.add_argument('-o','--out', dest = 'out', help = 'Output list file')
    args = parser.parse_args()
    import os
    args._in = os.path.realpath(args._in)
    if args.out != None: args.out = os.path.realpath(args.out)
    args.pheno.sort()
    
    from _utils import cmdhistory, logger
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()