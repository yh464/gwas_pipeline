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

def parse_greml_h2_log(file):
    '''
    Parses GREML H2 logs
    input: file name *.greml.hsq
    output: h2 and se
    '''
    try:
        tmp = open(file).read().splitlines()[4].split()
        h2 = float(tmp[-2])
        se = float(tmp[-1])
    except: h2 = np.nan; se = np.nan
    return h2, se

def parse_rg_log(file, full = False):
    '''
    Parses LDSC rg log for only one pair of phenotypes
    input: file name *.rg.log
    output: data frame with following columns: group1, pheno1, group2, pheno2, rg, se, p
    full output (specify full = True): in addition to above output, z scores;
        h2_obs, h2_int, gcov_int and their SE
    '''
    def floatna(x):
        try: return float(x)
        except: return np.nan
    hdr = ['group1','pheno1','group2','pheno2','rg','se','z','p','h2_obs',
           'h2_obs_se','h2_int','h2_int_se','gcov_int','gcov_int_se']
    all_stats = []
    tmp = open(file)
    skip = True
    line = 'placeholder'
    while len(line) > 0:
        line = tmp.readline()
        if line.find('gcov_int_se') > -1: skip = False; continue
        if line.find('Analysis finished') > -1: skip = True; continue
        if skip: continue
        tmp_stats = line.replace('\n','').split()
        if len(tmp_stats) == 0: continue
        group1 = os.path.basename(os.path.dirname(tmp_stats[0]))
        pheno1 = os.path.basename(tmp_stats[0]).replace('.sumstats','').replace('.gz','')
        group2 = os.path.basename(os.path.dirname(tmp_stats[1]))
        pheno2 = os.path.basename(tmp_stats[1]).replace('.sumstats','').replace('.gz','')
        tmp_stats = [group1, pheno1, group2, pheno2] + [floatna(x) for x in tmp_stats[2:]]
        all_stats.append(tmp_stats)
    
    if len(all_stats) == 0: return pd.DataFrame(data = [], index = [], columns = hdr)
    all_stats = pd.DataFrame(data = all_stats, columns = hdr)
    all_stats.loc[all_stats.rg > 1, 'rg'] = 1
    all_stats.loc[all_stats.rg < -1, 'rg'] = -1
    all_stats.loc[all_stats.se < 1e-20, 'se'] = 1e-20

    if full: return all_stats
    else: return all_stats[['group1','pheno1','group2','pheno2','rg','se','p']]