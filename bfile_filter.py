#!/usr/bin/env python3
'''
This programme filters the GRM to different thresholds
'''

import argparse

# argument input
parser = argparse.ArgumentParser(description=
  'This programme filters the GRM to different thresholds')
parser.add_argument('freq', help = 'Minor allele frequency threshold')
parser.add_argument('-b','--bdir', dest = 'bdir', help = 'BED file directory',
  default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Genetics/Genetic_data/Neuroimaging_samples/')
parser.add_argument('--plink', dest = 'plink', help = 'Location of plink2 executable',
  default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/plink2')
parser.add_argument('-o','--out',dest  = 'out', help = 'Output directory',
  default = '../params/')
parser.add_argument('-f','--force', dest = 'force', help = 'Force output',
  default = False, const = True, action = 'store_const')
args = parser.parse_args()

import os
args.out = os.path.realpath(args.out)

chrs = []
for tmp in range(1,23):
  chrs.append('ukbchr_v2_r2correct_v2_'+str(tmp))
chrs.append('ukbchr_v2_r2correct_X')
tmpdir = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/temp/'

if not os.path.isdir(tmpdir): os.system(f'mkdir -p {tmpdir}')

for c in chrs:
  if not os.path.isfile(f'{tmpdir}{c}.snplist'):
    out = os.system(f'{args.plink} --bfile {args.bdir}/{c} --max-maf {args.freq} --out {tmpdir}/{c}'+
              ' --write-snplist --rm-dup force-first')
    if out != 0: print(f'WARNING: snplist failed for chromosome {c}')
    print()

if os.path.isfile(f'{args.out}/{args.freq}.snplist'): 
  os.system(f'rm {args.out}/{args.freq}.snplist')
os.system(f'cat {tmpdir}/*.snplist >> {args.out}/{args.freq}.snplist')