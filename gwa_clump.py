#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2024-12-03

Clumps independent loci from a fastGWA format file

Requires following inputs: 
    GWAS summary statistics (single file)
'''

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i','--in', dest = '_in', help = 'Input directory')
parser.add_argument('--file', dest = 'file', help = 'Input file (fastGWA)')
parser.add_argument('-b','--bfile', dest = 'bfile', help = 'BED file list',
  default = '../params/bed_files_ukb.txt')
parser.add_argument('--plink', dest = 'plink', help = 'Path to PLINK *1.9* executable', 
  default = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Genetics/plink')
parser.add_argument('-o','--out', dest = 'out', help = 'Output directory')     # defaults to input dir
parser.add_argument('-p',help = 'p-value threshold',
  default = 5e-8, type = float) # or 3.1076e-11, or 5e-6; 3.1076e-11 is derived from matrix decomposition
parser.add_argument('-f','--force', dest = 'force', help = 'Force overwrite',
  default = False, action = 'store_true')
args = parser.parse_args()

import os
for arg in ['_in','out','bfile']:
    exec(f'args.{arg} = os.path.realpath(args.{arg})')
if type(args.out) == type(None): args.out = args._in

from _utils import logger
logger.splash(args)

import time
import pandas as pd
import numpy as np

tic = time.perf_counter()
idx = 0
blist = np.loadtxt(args.bfile,dtype = 'U')
prefix = '.'.join(args.file.split('.')[:-1])
out = f'{args.out}/{prefix}_{args.p:.0e}.clumped'

tmpdir = f'/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/temp/clump_cache/{os.path.basename(args._in)}_{args.p:.0e}'
logdir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/logs/gwa_clump/'
if not os.path.isdir(tmpdir): os.system(f'mkdir -p {tmpdir}')
if not os.path.isdir(logdir): os.mkdir(logdir)
log = open(f'{logdir}/{args.file}.log'.replace('.fastGWA',''),'w')     # log file

os.chdir(args.out)                                                             # we do not need the input dir

df = pd.read_table(f'{args._in}/{args.file}', sep = '\t')
sf = df.P.values < args.p                                                      # sig filter, must be determined by matrix decomposition
if sf.sum() == 0:
  print(f'File {args.file} contains no significant SNP, skipping', file = log)
  toc = time.perf_counter() - tic
  print(f'Total time = {toc:.3f}.', file = log)

df_sig = df.loc[sf,:].sort_values(by = 'P')
df_sig.to_csv(out.replace('clumped','siglist'), sep = '\t',index = False)      # export top few snps

tmp_flist = []                                                                 # list of temp files
chrs = df_sig['CHR'].unique()
idx = 0

out_df = []
for c in chrs:                                                                 # every chromosome that is sig, sorted 
  idx += 1
  # separates files by chromosome
  df_tmp = df.loc[df.CHR == c, :].sort_values(by = 'P')
  tmpgwa = f'{tmpdir}/{args.file}_chr{c}.fastGWA'
  tmpsnp = f'{tmpdir}/{args.file}_chr{c}.snplist'
  df_tmp.to_csv(tmpgwa, index = False, sep = '\t')
  df_tmp['SNP'].to_csv(tmpsnp, index = False, header = False)
  bf = blist[c-1]                                                              # c ranges 1-23
  tmpout = f'{tmpdir}/{args.file}_chr{c}'
  tmp_flist.append(tmpgwa)
  tmp_flist.append(tmpsnp)
  tmp_flist.append(f'{tmpout}.hh')
  tmp_flist.append(f'{tmpout}.clumped')
  
  if not os.path.isfile(f'{tmpout}.clumped') or args.force:
      # clumps by chromosome
      os.system(f'{args.plink} --noweb --bfile {bf} --clump {tmpgwa} '+
        f'--clump-field P --clump-p1 {args.p} --clump-p2 1 --clump-r2 0.1 '+       # p1 must be determined by matrix decomposition
        f'--clump-kb 1000 --extract {tmpsnp} --out {tmpout}')
  
  out_df.append(pd.read_table(f'{tmpout}.clumped', sep = '\s+'))
  toc = time.perf_counter() - tic
  print(f'Finished clumping chromosome {c}, {idx}/{len(chrs)} time = {toc:.3f}.', file = log)

out_df = pd.concat(out_df, axis = 0)
out_df.to_csv(out, sep = '\t', index = False)
# # clears temp
# for x in tmp_flist:
#   try: os.remove(x)
#   except: pass