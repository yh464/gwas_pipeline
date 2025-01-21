#!/usr/bin/env python3
'''
This script constitutes the LDSC heritability pipeline for any incoming fastGWA file
'''

import argparse
parser = argparse.ArgumentParser(description = 
  'This script computes heritability analysis for any incoming fastGWA file')
parser.add_argument('--ldsc', dest = 'ldsc', help = 'LDSC executable directory',
  default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/ldsc/')
parser.add_argument('-i', '--in', dest = '_in', help = 'input fastGWA file')
parser.add_argument('-o','--out', dest = 'out', help = 'output directory (ABSOLUTE)')
parser.add_argument('-f','--force',dest = 'force', help = 'force output',
  default = False, action = 'store_true')
args = parser.parse_args()

import os

args.out = os.path.realpath(args.out)
prefix = os.path.basename(args._in).replace('.fastGWA', '')

scripts_path = os.path.realpath(__file__)
scripts_path = os.path.dirname(scripts_path)
if args.force or (not os.path.isfile(f'{args.out}/{prefix}.sumstats')):
    # this command uses python2 so a separate script for ldsc  
    os.system(f'bash {scripts_path}/ldsc_master.sh munge_sumstats.py --sumstats {args._in} '+ \
          # f'--merge-alleles {args.ldsc}/ukb_snp_info.txt '+
          f'--out {args.out}/{prefix}')

# QC h2 log
h2log = f'{args.out}/{prefix}.h2.log'
if os.path.isfile(h2log):
    tmp = open(h2log).read().splitlines()[-7].replace('(','').replace(')','').split()
    try: h2 = float(tmp[-2])
    except: os.remove(h2log)

if (not os.path.isfile(f'{args.out}/{prefix}.h2.log')) or args.force:
  os.system(f'bash {scripts_path}/ldsc_master.sh ldsc.py '+
    f'--ref-ld-chr {args.ldsc}/baseline/ --w-ld-chr {args.ldsc}/baseline/ '+
    f'--h2 {args.out}/{prefix}.sumstats '+
    f'--out {args.out}/{prefix}.h2')