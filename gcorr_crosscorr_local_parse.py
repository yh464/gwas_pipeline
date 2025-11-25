#!/usr/env/bin python
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1.0: 2023-07-15
Version 1.1: 2024-12-09
Version 1.2: 2025-02-11

Creates a summary table for genetic correlations for local phenotypes

Upstream workflow: gcorr_crosscorr_local_batch.py
Input format: LONG-format table, group1 (*local) pheno1 group2 pheno2 rg se p q
'''

import argparse
parser = argparse.ArgumentParser(
  description = 'this script prepares the psychiatric rg for r-ggseg plotting')
parser.add_argument('-p1', nargs = '*', help = 'local phenotypes',
  default=['deg_local','degi_local','degc_local','clu_local','eff_local','mpl_local'])
parser.add_argument('-p2', nargs = '*', help = 'correlates phenotype groups',
  default = ['disorders'])
parser.add_argument('-i','--in', dest = '_in', help = 'input directory',
  default='../local_corr/')
parser.add_argument('--sumstats', help = 'sumstats directory to be scanned for file names',
  default = '../gcorr/ldsc_sumstats/')
parser.add_argument('-c','--corresponding', dest = 'corresponding', 
  help = 'match local phenotype to corresponding global phenotype',
  default = False, action = 'store_true')
parser.add_argument('-o','--out', dest = 'out', help = 'output directory')
# always overwrites
args = parser.parse_args()

from ._utils import logger, path
logger.splash(args)
norm = path.normaliser()

import os
import numpy as np
import pandas as pd
import scipy.stats as sts
from fnmatch import fnmatch
args._in = os.path.realpath(args._in)
args.out = args._in if type(args.out) == type(None) else os.path.realpath(args.out)

# scan for file names
os.chdir(args.sumstats)
prefix_1 = []; pheno_1 = []
prefix_2 = []; pheno_2 = []
for p in args.p1:
    for x in os.listdir(p):
        if fnmatch(x,'*.sumstats'):
            prefix_1.append(x.replace('.sumstats','')); pheno_1.append(p)
for p in args.p2:
    for x in os.listdir(p):
        if fnmatch(x,'*.sumstats'):
            prefix_2.append(x.replace('.sumstats','')); pheno_2.append(p)

os.chdir(args._in)
summary = []
for g2, p2 in zip(pheno_2, prefix_2):
    for g1, p1 in zip(pheno_1, prefix_1):
        if args.corresponding:
            if p2.find(g1.replace('local','')) == -1: continue
            if fnmatch(p2, '*_l_*') or fnmatch(p2, '*_r_*') or \
                fnmatch(p2, '*_l') or fnmatch(p2, '*_r'): 
                continue
        fname = f'{args._in}/{g1}/{g2}/{g1}_{p1}.{g2}_{p2}.rg.log'
        if not os.path.isfile(fname):
            print(f'{fname} does not exist')
            continue
        tmp = open(fname)
        tmp_stats = tmp.read().splitlines()
        tmp_stats = tmp_stats[-4].split()
        while tmp_stats.count('') > 0:
          tmp_stats.remove('')
        try: 
            rg = float(tmp_stats[2])
            if rg > 1: rg = 1
            if rg < -1: rg = -1
        except: 
          rg = np.nan
          print(f'{fname} shows NA correlation!')
        try: se = max((float(tmp_stats[3]),10**-20))
        except: se = np.nan

        p = 1-sts.chi2.cdf((rg/se)**2, df = 1) # p value
        
        summary.append(pd.DataFrame(dict(group1 = g1, pheno1 = p1,
          group2 = g2, pheno2 = p2, rg = rg, se = se, p = [p])))

summary = pd.concat(summary)

for g2 in args.p2:
    if args.corresponding:
        temp = summary.loc[summary.group2 == g2,:]
        norm.normalise(temp).to_csv(f'{args.out}/gcorr_local_{g2}.long.txt', sep = '\t', 
                    index = False, header = True)
        outdf = temp.pivot(columns = 'group1',index = 'pheno1', values = 'rg')
        outdf.columns.name = None
        outdf.index.name = 'label'
        norm.normalise(outdf).to_csv(f'{args.out}/gcorr_local_{g2}.txt', sep = '\t', index = True, header = True)
    else:
        for z in summary.pheno2.unique():
            outdf = summary.loc[summary.pheno2 == z,:]
            norm.normalise(outdf).to_csv(f'{args.out}/gcorr_local_{z}.long.txt', sep = '\t', 
                         index = False, header = True)
            outdf = outdf.pivot(columns = 'group1', index = 'pheno1', values = 'rg')
            outdf.columns.name = None
            outdf.index.name = 'label'
            norm.normalise(outdf).to_csv(f'{args.out}/gcorr_local_{z}.txt', 
                                         sep = '\t', index = True, header = True)
