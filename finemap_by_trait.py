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

import argparse
parser = argparse.ArgumentParser(
  description = 'This programme constitutes the finemap pipeline for a phenotype')
parser.add_argument('pheno', help = 'Phenotype')
parser.add_argument('-i','--in', dest = '_in', help = 'Input GWA p-statistic file')
parser.add_argument('-d','--dir', dest = 'dir', help = 'Directory containing all GWAS summary stats',
  default = '../gwa/')
parser.add_argument('-c','--clump', dest = 'clump', help = 'Directory containing all clump outputs',
  default = '../clump/')
parser.add_argument('-o', '--out', dest = 'out', help = 'output directory',
  default = '../finemap/')
parser.add_argument('-b', '--bfile', dest = 'bfile', help = 'directory of bed binaries',
  default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/bed/')
parser.add_argument('--polyfun', help = 'directory of POLYFUN tool',
  default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/polyfun/')
parser.add_argument('-p', dest = 'p', help = 'p-value', type = float, default = 3.1076e-11)
parser.add_argument('-f','--force',dest = 'force', help = 'force output',
  default = False, action = 'store_true')
args = parser.parse_args()

import os
for arg in ['dir','clump','out','bfile']:
    setattr(args, arg, os.path.realpath(getattr(args, arg)))

from _utils import logger
logger.splash(args)

if not os.path.isdir(args.out): os.mkdir(args.out)
os.chdir(args.out)
if not os.path.isdir(args.pheno): os.mkdir(args.pheno)
os.chdir(args.pheno)
stats_dir = f'{args.out}/{args.pheno}/polyfun_stats'
if not os.path.isdir(stats_dir): os.system(f'mkdir -p {stats_dir}')

os.chdir(f'{args.dir}/{args.pheno}')
prefix = args._in.replace('.fastGWA','')

# need to specify the environment of python - use absolute path of interpreter
o = 0
if (not os.path.isfile(f'{stats_dir}/{prefix}.polyfun.stats')) or args.force:
  o = os.system('python '+
            f'{args.polyfun}/munge_polyfun_sumstats.py '+
            f'--sumstats {args.dir}/{args.pheno}/{args._in} --out {stats_dir}/{prefix}.polyfun.stats')
  if o != 0: raise Exception(f'ERROR for {args._in} at step 1: munge_polyfun_sumstats')

o = 0
if (not os.path.isfile(f'{stats_dir}/{prefix}.snpvar')) or args.force:
  o = os.system('python '+
            f'{args.polyfun}/extract_snpvar.py --allow-missing '+
            f'--sumstats {stats_dir}/{prefix}.polyfun.stats --out {stats_dir}/{prefix}.snpvar')
  if o != 0: raise Exception(f'ERROR for {args._in} at step 2: extract snpvar')

import pandas as pd
# identify SNPs that need to be clumped
overlaps = pd.read_table(f'{args.clump}/{args.pheno}_{args.p:.0e}_overlaps.txt', index_col ='label').T
prefix = args._in.replace('.fastGWA','')
overlaps = overlaps[prefix]
overlaps = overlaps[overlaps > 0]
snps = overlaps.index

df = pd.read_csv(f'{args.dir}/{args.pheno}/{args._in}', sep = '\s+', 
                 usecols = ['SNP','CHR','POS','N'])
n = df['N'].max()

summary = []
print('Following SNPs are being fine-mapped')
print(snps.to_numpy())
for snp in snps:
    pos = df.loc[df.SNP == snp, 'POS'].iloc[0]
    start = int(min(pos-5*10**5,1))
    stop = start + 10**6
    c = df.loc[df.SNP == snp, 'CHR'].iloc[0]
    if c < 23: geno = f'{args.bfile}/chr{c:.0f}'
    else: geno = f'{args.bfile}/chrX'
    
    # finemap +/- 0.5 mb for causal genes
    o = 0
    if (not os.path.isfile(f'{stats_dir}/{prefix}_chr{c}_{start}_{stop}.csv')) or args.force:
      scripts_path = os.path.realpath(__file__)
      scripts_path = os.path.dirname(scripts_path)
      cmd = f'python {args.polyfun}/finemapper.py '+ \
        f'--method susie --n {n:.0f} --sumstats {stats_dir}/{prefix}.snpvar --chr {c:.0f} ' + \
        f'--start {start:.0f} --end {stop:.0f} --geno {geno} --out {stats_dir}/{prefix}_chr{c}_{start}_{stop}.txt '+ \
        '--max-num-causal 5 --allow-swapped-indel-alleles'
      o = os.system(cmd)
      
      # o = os.system(f'bash {scripts_path}/'+
      #   f'finemap_by_seg.sh {stats_dir}/{prefix}.snpvar {c:.0f} {start:.0f} {stop:.0f} {geno} '+
      #   f'{stats_dir}/{prefix}_chr{c}_{start}_{stop}.txt')
      if o != 0: 
          print(cmd)
          continue
          # raise Exception(f'ERROR for {args._in} at step 3: finemap {c}/{start}/{stop}')
    
    summary.append(pd.read_table(f'{stats_dir}/{prefix}_chr{c}_{start}_{stop}.txt'))

summary = pd.concat(summary)
summary.to_csv(f'{args.out}/{args.pheno}/{prefix}.finemap.summary', sep = '\t',
               index = False)