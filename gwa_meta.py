#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2024-11-29
Version 2: 2024-12-03

Conducts and formats GWAS meta-analysis

Requires following inputs:
    GWAS summary stats
Changelog:
    uses PLINK to estimate beta, SE and Z instead of METAL
'''

def main(args):
    import os
    import pandas as pd
    import scipy.stats as sts
    
    # output format
    out_prefix = args.out.replace('.fastGWA','')
    
    # METAL: heterogeneity test, AF1/SE, 
    metal_out = f'{out_prefix}.metal'
    if not os.path.isfile(metal_out) or args.force:
        tmpdir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/temp/metal_scripts/'
        if not os.path.isdir(tmpdir): os.system(f'mkdir -p {tmpdir}')
        
        # temp metal script name
        finfo = []
        for x in args._in:
            finfo.append(os.path.basename(os.path.dirname(x)))
        finfo.append(os.path.basename(x))
        mts = f'{tmpdir}/'+ '_'.join(finfo) + '.metal'
        mts = open(mts, 'w')
        
        # standard metal settings
        print('marker SNP', file = mts)
        print('allele A1 A2', file = mts)
        print('effect BETA', file = mts)
        print('pvalue P', file = mts)
        print('freq AF1', file = mts)
        print('averagefreq on', file = mts)
        
        # specify output
        print(f'out {out_prefix}_ .fastGWA', file = mts)
        
        # specify input
        for x in args._in:
            print(f'process {x}', file = mts)
            
        print('analyze heterogeneity', file = mts)
        print('quit', file = mts)
        mts.close() # IMPORTANT!
        
        os.system(f'{args.metal} {tmpdir}/'+ '_'.join(finfo) + '.metal')
        os.system(f'mv {out_prefix}_1.fastGWA {metal_out}')
        os.remove(f'{out_prefix}_1.fastGWA.info')
    
    # PLINK: BETA 
    plink_out = f'{out_prefix}.meta'
    if not os.path.isfile(plink_out) or args.force:
        os.system(f'{args.plink} --meta-analysis '+ ' '.join(args._in) +
            f' + qt --meta-analysis-bp-field POS --out {out_prefix}')
    
    # combine results to produce output file
    if not os.path.isfile(args.out) or args.force:
        # process METAL outputs
        metal_df = pd.read_table(metal_out)
        # format columns
        metal_df.columns = ['SNP','A1x','A2','AF1','AF1SE','N','Z','P','DIR',
                            'HET_I2','HET_CHI2','HET_DF','HET_P']
        metal_df['A1x'] = metal_df['A1x'].str.upper()
        metal_df = metal_df[['SNP','A1x','AF1','AF1SE','N','DIR',
                             'HET_I2','HET_CHI2','HET_DF','HET_P']]
        print(metal_df.head())
        # estimate BETA and SE from PLINK output
        plink_df = pd.read_table(plink_out, sep = '\s+')
        plink_df = plink_df[['CHR','SNP','BP','A1','A2','P','BETA']]
        plink_df.columns = ['CHR','SNP','POS','A1','A2','P','BETA']
        plink_df['SE'] = abs(plink_df['BETA'] / sts.norm.ppf(plink_df['P']/2))
        plink_df['Z'] = plink_df['BETA']/plink_df['SE']
        plink_df.dropna(inplace=True)
        print(plink_df.head())
        # harmonise allele frequencies
        df = pd.merge(plink_df, metal_df, on = ['SNP'])
        df.loc[df['A1x'] != df['A1'], 'AF1'] = 1 - df.loc[df['A1x'] != df['A1'], 'AF1']
        df.drop('A1x', axis = 1, inplace = True)
        df = df[['CHR','SNP','POS','A1','A2','AF1','AF1SE','N','BETA','SE','Z','P',
                 'DIR','HET_I2','HET_CHI2','HET_DF','HET_P']].dropna()
        # sort by p values
        df = df.sort_values(by = ['CHR','POS'], ascending = True)
        df.to_csv(args.out, sep = '\t', index = False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'This programme creates genetic correlation matrices for global phenotypes')
    parser.add_argument('-i','--in', dest = '_in', help = 'input files in fastGWA format',
      nargs = '*')
    parser.add_argument('--metal', help = 'METAL executable',
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/metal') # intended to be absolute
    parser.add_argument('--plink', help = 'PLINK 1.9 executable',
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/plink') # intended to be absolute
    parser.add_argument('-o','--out', dest = 'out', help = 'output file name')
    parser.add_argument('-f','--force',dest = 'force', help = 'force overwrite',
      default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['out','metal','plink']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
    from _utils import cmdhistory, logger
    logger.splash(args)
    main(args)