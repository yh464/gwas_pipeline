#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-03-24

Parser functions for genetic logs (IMPORTANT)
This script is specific to the GWAS pipeline and will not be migrated to _utils
'''

import os
from fnmatch import fnmatch
import pandas as pd
import numpy as np
import warnings

def find_gwas(*pheno, 
              dirname = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/gwa', 
              ext = 'fastGWA'):
    '''
    Data structure: {dirname}/{pheno[0]}/*.{ext}
    pheno: phenotype groups
    dirname: directory of all GWAS sumstats
    ext: extension, usually fastGWA
    '''
    out = []
    for p in pheno:
        for x in os.listdir(f'{dirname}/{p}'):
            if not fnmatch(x, f'*.{ext}') or fnmatch (x,f'*_X.{ext}'): continue
            out.append((p,x.replace(f'.{ext}','')))
    return out

def find_clump(dirname, prefix, pval):
    '''
    Find PLINK clump files for a specific trait
    Quality controls to find strictest p-value threshold with >5 SNP
    dirname: Directory to look for clumps
    prefix: name of phenotype
    pval: p-value
    '''
    if os.path.isfile(f'{dirname}/{prefix}_{pval:.0e}.clumped'):
        # min 5 SNPs
        if len(open(f'{dirname}/{prefix}_{pval:.0e}.clumped').read().splitlines()) > 5:
            return f'{dirname}/{prefix}_{pval:.0e}.clumped', pval
    # identify clump file with lowest p-value with >=5 SNPs
    flist = [] 
    for y in os.listdir(dirname):
        if fnmatch(y,f'{prefix}_?e-??.clumped'): 
            if len(open(f'{dirname}/{y}').read().splitlines()) > 5:
                flist.append(y)
    if len(flist) > 0:
        plist = [float(z[-13:-8]) for z in flist]
        return f'{dirname}/{prefix}_{min(plist):.0e}.clumped', min(plist)
    
    for y in os.listdir(dirname):
        if fnmatch(y,f'{prefix}_?e-??.clumped'): 
             flist.append(y)
    if len(flist) > 0:
        plist = [float(z[-13:-8]) for z in flist]
        warnings.warn(f'{prefix} has <5 SNPs')
        return f'{dirname}/{prefix}_{max(plist):.0e}.clumped', max(plist)
    raise FileNotFoundError(f'No clump found for {prefix}')
    
def parse_h2_log(file):
    '''
    Parses LDSC H2 logs
    input: file name *.h2.log
    output: h2 and se
    '''
    try:
        for line in open(file).read().splitlines():
            if fnmatch(line, 'Total Observed scale h2*'): break
        l = line.split()
        h2 = float(l[-2])
        se = float(l[-1].replace('(','').replace(')',''))
    except:
        h2 = np.nan; se = np.nan
    return h2, se

def parse_rg_log(file):
    '''
    Parses LDSC H2 logs
    input: file name *.rg.log
    output: rg and se
    '''
    tmp = open(file)
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
      print(f'{file} shows NA correlation!')
    try: se = max((float(tmp_stats[3]),10**-20))
    except: se = np.nan

    return rg, se