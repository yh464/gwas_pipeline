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
parser.add_argument('-b', '--bfile', dest = 'bfile', help = 'bed binary',
  default = '../params/bed/')
parser.add_argument('-p', dest = 'p', help = 'p-value', type = float, default = 3.1076e-11)
parser.add_argument('-f','--force',dest = 'force', help = 'force output',
  default = False, action = 'store_true')
args = parser.parse_args()

import os
for arg in ['dir','out','bfile']:
    exec(f'args.{arg} = os.path.realpath(args.{arg})')

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
  o = os.system('/home/yh464/.conda/envs/gentoolspy/bin/python '+
            '/rds/user/yh464/hpc-work/conda/polyfun/munge_polyfun_sumstats.py '+
            f'--sumstats {args.dir}/{args.pheno}/{args._in} --out {stats_dir}/{prefix}.polyfun.stats')
  if o != 0: raise Exception(f'ERROR for {args._in} at step 1: munge_polyfun_sumstats')

o = 0
if (not os.path.isfile(f'{stats_dir}/{prefix}.snpvar')) or args.force:
  o = os.system('/home/yh464/.conda/envs/gentoolspy/bin/python '+
            '/rds/user/yh464/hpc-work/conda/polyfun/extract_snpvar.py --allow-missing '+
            f'--sumstats {stats_dir}/{prefix}.polyfun.stats --out {stats_dir}/{prefix}.snpvar')
  if o != 0: raise Exception(f'ERROR for {args._in} at step 2: extract snpvar')

import pandas as pd
df = pd.read_csv(f'{args.dir}/{args.pheno}/{args._in}', sep = '\s+').sort_values(
  by = 'P', ignore_index = True)

# repeat the fine-mapping for all SNPs with p<3.1076e-11
idx = 0
while True:
  p = df['P'][idx]
  if p > args.p: break
  pos = df['POS'][idx]
  start = pos-5*10**5
  stop = pos+5*10**5
  c = df['CHR'][idx]
  if c < 23: geno = f'{args.bfile}/chr{c:.0f}'
  else: geno = f'{args.bfile}/chrX'
  
  # finemap +/- 0.5 mb for causal genes
  o = 0
  if (not os.path.isfile(f'{stats_dir}/{prefix}_chr{c}_{start}_{stop}.csv')) or args.force:
    scripts_path = os.path.realpath(__file__)
    scripts_path = os.path.dirname(scripts_path)
    o = os.system(f'bash {scripts_path}/'+
      f'finemap_by_seg.sh {stats_dir}/{prefix}.snpvar {c:.0f} {start:.0f} {stop:.0f} {geno} '+
      f'{stats_dir}/{prefix}_chr{c}_{start}_{stop}.csv')
    if o != 0: raise Exception(f'ERROR for {args._in} at step 3: finemap {c}/{start}/{stop}')
  
  # concatenate results for reference
  if idx == 0:
    os.system(f'cat {stats_dir}/{prefix}_chr{c}_{start}_{stop}.csv > '+
              f'{args.out}/{args.pheno}/{prefix}.finemap.summary')
  else:
    os.system(f'tail -n +2 {stats_dir}/{prefix}_chr{c}_{start}_{stop}.csv >> '+
              f'{args.out}/{args.pheno}/{prefix}.finemap.summary')
  idx += 1