#!/usr/env/bin python
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1.0: 2023-07-15
Version 1.1: 2024-12-09

Creates a summary table for genetic correlations for local phenotypes

Upstream workflow: gcorr_crosscorr_plot.py
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
parser.add_argument('-c','--corresponding', dest = 'corresponding', 
  help = 'match local phenotype to corresponding global phenotype',
  default = False, action = 'store_true')
parser.add_argument('-o','--out', dest = 'out', help = 'output directory')
# always overwrites
args = parser.parse_args()

from _utils import logger
logger.splash(args)

import os
args._in = os.path.realpath(args._in)
args.out = args._in if type(args.out) == type(None) else os.path.realpath(args.out)
os.chdir(args._in)

import pandas as pd
from fnmatch import fnmatch # to detect if we are correlating local v global
for y in args.p2:    
    dflist = []
    for x in args.p1:
        tmp = pd.read_table(f'{x}/crosscorr_{x}.{y}.txt')
        if args.corresponding: # filter input df
            # identify matching pheno2 labels
            t1 = []
            for z in tmp.pheno2.unique():
                if fnmatch(z,'*'+x.replace('local','')+'*'): t1.append(z)
            tmp = tmp.loc[tmp.pheno2.isin(t1),:]
        dflist.append(tmp)
    summary = pd.concat(dflist)
    print(summary.head())
    
    if args.corresponding:
        outdf = summary.pivot(columns = 'group1',index = 'pheno1', values = 'rg')
        outdf.columns.name = None
        outdf.index.name = 'label'
        outdf.to_csv(f'{args.out}/{y}.txt', sep = '\t', index = True, header = True)
    else:
        for z in summary.pheno2.unique():
            outdf = summary.loc[summary.pheno2 == z,:].pivot(
                columns = 'group1', index = 'pheno1', values = 'rg')
            outdf.columns.name = None
            outdf.index.name = 'label'
            outdf.to_csv(f'{args.out}/{z}.txt', sep = '\t', index = True, header = True)
