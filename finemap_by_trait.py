#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1.0: 2023-07-22
Version 1.1: 2024-12-04

This is a script for single-trait fine-mapping using polyfun-susie

Preceding workflow:
    gwa_batch.py
    gwa_clump_batch.py
    gwa_clump_parse.py
Required input:
    GWAS summary statistics (single file)
    
Changelog:
    Now requires clumping output to reduce unnecessary computation
'''

def main(args):
    if os.path.isfile(args.out) and (not args.force): return

    stats_dir = f'{os.path.dirname(args.out)}/polyfun_stats'
    if not os.path.isdir(stats_dir): os.system(f'mkdir -p {stats_dir}')
    prefix = os.path.basename(args._in).replace('.fastGWA','')

    # need to specify the environment of python - use absolute path of interpreter
    o = 0
    if (not os.path.isfile(f'{stats_dir}/{prefix}.polyfun.stats')) or args.force:
      o = os.system('python '+
                f'{args.polyfun}/munge_polyfun_sumstats.py '+
                f'--sumstats {args._in} --out {stats_dir}/{prefix}.polyfun.stats')
      if o != 0: raise Exception(f'ERROR for {args._in} at step 1: munge_polyfun_sumstats')

    o = 0
    if (not os.path.isfile(f'{stats_dir}/{prefix}.snpvar')) or args.force:
      o = os.system('python '+
                f'{args.polyfun}/extract_snpvar.py --allow-missing '+
                f'--sumstats {stats_dir}/{prefix}.polyfun.stats --out {stats_dir}/{prefix}.snpvar')
      if o != 0: raise Exception(f'ERROR for {args._in} at step 2: extract snpvar')

    import pandas as pd
    from ._plugins.logparser import parse_clump_file, overlap_clumps
    from ._utils.path import find_bed
    # identify SNPs that need to be clumped
    clumps = parse_clump_file(args.clump).sort_values(['CHR','POS']).reset_index(drop=True)
    clumps['group'] = prefix; clumps['pheno'] = os.path.basename(os.path.dirname(args._in))
    overlaps,_ = overlap_clumps(clumps, dist = 5e6)
    snps = overlaps.SNP

    # identify overlaps
    print('Following SNPs are being fine-mapped')
    print(snps.to_numpy())

    # sample size
    df = pd.read_csv(args._in, sep = '\\s+', usecols = ['N'])
    n = df['N'].max()

    # find PLINK binaries
    bed = find_bed(args.bfile, x = True)

    summary = []
    for i in range(overlaps.shape[0]):
        start = int(overlaps.START.iloc[i])
        stop = int(overlaps.STOP.iloc[i])
        c = int(overlaps.CHR.iloc[i])
        geno = bed[c-1]
        
        # finemap +/- 0.5 mb for causal variants
        o = 0
        if (not os.path.isfile(f'{stats_dir}/{prefix}_chr{c}_{start}_{stop}.txt')) or args.force:
          cmd = f'python {args.polyfun}/finemapper.py '+ \
            f'--method susie --n {n:.0f} --sumstats {stats_dir}/{prefix}.snpvar --chr {c} ' + \
            f'--start {start:.0f} --end {stop:.0f} --geno {geno} --out {stats_dir}/{prefix}_chr{c}_{start}_{stop}.txt '+ \
            '--max-num-causal 5 --allow-swapped-indel-alleles --allow-missing'
          o = os.system(cmd)

          if o != 0: 
              print(cmd)
              continue
        summary.append(pd.read_table(f'{stats_dir}/{prefix}_chr{c}_{start}_{stop}.txt'))

    summary = pd.concat(summary)
    summary.to_csv(args.out, sep = '\t', index = False)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
      description = 'This programme constitutes the finemap pipeline for a phenotype')
    parser.add_argument('-i','--in', dest = '_in', help = 'Input GWA p-statistic file', required = True)
    parser.add_argument('-c','--clump', dest = 'clump', help = 'Input clumping output file', required = True)
    parser.add_argument('-o', '--out', dest = 'out', help = 'output file name', required = True)
    parser.add_argument('-b', '--bfile', dest = 'bfile', help = 'directory of bed binaries',
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/bed/')
    parser.add_argument('--polyfun', help = 'directory of POLYFUN tool',
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/polyfun/')
    parser.add_argument('-p', dest = 'p', help = 'p-value', type = float, default = 5e-8)
    parser.add_argument('-f','--force',dest = 'force', help = 'force overwrite',
      default = False, action = 'store_true')
    args = parser.parse_args()

    import os
    for arg in ['_in','clump','out','bfile']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))

    from ._utils import logger
    logger.splash(args)

    main(args)