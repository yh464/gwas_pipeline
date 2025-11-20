#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-03-24

Parser functions for genetic logs (IMPORTANT)
This script is specific to the GWAS pipeline and will not be migrated to _utils
'''

import os, warnings
from fnmatch import fnmatch
import pandas as pd
import numpy as np
import scipy.stats as sts
    
def parse_h2_log(file, full = False):
    '''
    Parses LDSC H2 logs
    input: file name *.h2.log
    output: h2 and se
    '''
    try:
        for line in open(file).read().splitlines():
            if fnmatch(line, 'Total Observed scale h2*'): 
                l = line.split()
                h2 = float(l[-2])
                se = float(l[-1].replace('(','').replace(')',''))
            if fnmatch(line, 'Intercept*'):
                l = line.split()
                gcov_int = float(l[-2]); gcov_int_se = float(l[-1].replace('(','').replace(')',''))
        h2; se # important in case the total observed scale h2 is not read from log
    except:
        h2 = np.nan; se = np.nan; gcov_int = np.nan; gcov_int_se = np.nan
    if not full: return h2, se
    df = pd.DataFrame(dict(
        group1 = [os.path.basename(os.path.dirname(file))],
        pheno1 = os.path.basename(file).replace('.h2.log','').replace('.gz',''),
        group2 = os.path.basename(os.path.dirname(file)),
        pheno2 = os.path.basename(file).replace('.h2.log','').replace('.gz',''),
        rg = h2, se = se, p = 1 - sts.chi2.cdf(h2**2/se**2, df = 1),
        h2_obs = h2, h2_obs_se = se, h2_int = gcov_int, h2_int_se = gcov_int_se,
        gcov_int = gcov_int, gcov_int_se = gcov_int_se, z = h2/se
    ))
    return df

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

def crosscorr_parse(gwa1, gwa2 = [], 
        logdir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/gcorr/rglog',
        h2dir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/gcorr/ldsc_sumstats',
        exclude = [], full = False):
    '''
    gwa1 and gwa2 are lists of (group, pheno_list) tuples or (group, pheno) tuples, compatible with long/short
    leave gwa2 blank to estimate auto-correlations of gwa1
    '''
    summary = []
    
    from _utils.path import pair_gwas
    pairwise = pair_gwas(gwa1, gwa2)
    
    for g1, p1s, g2, p2s in pairwise:
        if g1 > g2: 
            g1, p1s, g2, p2s = g2, p2s, g1, p1s
            flip = True
        else: flip = False
        if isinstance(p1s, str): p1s = [p1s]
        if isinstance(p2s, str): p2s = [p2s]
        for p1 in p1s:
            if g1 == g2 and h2dir != None: # heritability
                fname = f'{h2dir}/{g1}/{p1}.h2.log'
                if not full:
                    h2, se = parse_h2_log(fname)
                    rg = pd.DataFrame(dict(group1 = [g1], pheno1 = p1, group2 = g1, 
                        pheno2 = p1, rg = h2, se = se, 
                        p = 1-sts.chi2.cdf(h2**2/se**2, df = 1),fixed_int = False))
                else: rg = parse_h2_log(fname, full = True)
                rg['fixed_int'] = False
                if flip: rg.iloc[:,[0,1,2,3]] = rg.iloc[:,[2,3,0,1]]
                summary.append(rg)
                
            fname = f'{logdir}/{g1}.{g2}/{g1}_{p1}.{g2}.rg.log'
            if not os.path.isfile(fname) and g1 != g2: 
                warnings.warn(f'No gene correlation found for {g1}/{p1} with {g2}\n'+
                    f'Try running:\n\n python gcorr_batch.py -p1 {g1} -p2 {g2}\n')
                continue
            elif not os.path.isfile(fname) and g1 == g2: continue
            rg = parse_rg_log(fname, full = full)
            rg['fixed_int'] = False
            rg = rg.loc[rg.pheno2.isin(p2s),:] # due to the file structure, some other traits may be present in the rg log
            if flip: rg.iloc[:,[0,1,2,3]] = rg.iloc[:,[2,3,0,1]]
            summary.append(rg)
            
            fname_noint = fname.replace('.rg.log','.noint.rg.log')
            if os.path.isfile(fname_noint):
                rg = parse_rg_log(fname_noint, full = full)
                rg['fixed_int'] = True
                if flip: rg.iloc[:,[0,1,2,3]] = rg.iloc[:,[2,3,0,1]]
                summary.append(rg)
            
    summary = pd.concat(summary) # creates a long format table
    summary.insert(loc = len(summary.columns), column = 'q', value = np.nan)
    for g1,p1s in gwa1: # FDR correction for each IDP, which are non-independent
        if isinstance(p1s,str): p1s = [p1s]
        for p1 in p1s:
            summary.loc[(summary.pheno1==p1) & (summary.group1==g1) & ~np.isnan(summary.p),'q'] \
                = sts.false_discovery_control(
            summary.loc[(summary.pheno1==p1)&(summary.group1==g1) & ~np.isnan(summary.p),'p'])
    if len(exclude) > 0:
        exclude_g = [x.split('/')[0] for x in exclude]
        exclude_p = [x.split('/')[1] for x in exclude]
        summary = summary.loc[~((summary.group1.isin(exclude_g) & summary.pheno1.isin(exclude_p)) |
            (summary.group2.isin(exclude_g) & summary.pheno2.isin(exclude_p))),:]
    return summary.drop_duplicates().reset_index(drop = True)

def parse_clump_file(file):
    return pd.read_table(file, sep = '\\s+').rename(columns = {'BP':'POS'})

class clump():
    # automatically identifies overlapping clumps from dataframe ROWS 
    # from parse_clump_file with added 'group' and 'pheno'
    def __init__(self, entry):
        self.chr = entry['CHR']
        self.start = entry['POS']; self.stop = entry['POS']
        self.sig_snps = {entry['SNP']}
        self.sentinel_snp = entry['SNP']
        self.pval = entry['P']
        self.phenotypes = {(entry['group'], entry['pheno'])}
        self.snps = set(entry['SP2'].replace('(1)','').split(','))

    # returns False if the new entry is independent from the clump
    def update(self, entry, dist:int = 0):
        if entry['CHR'] != self.chr: return False
        if not dist and entry['SNP'] not in self.snps: return False
        if dist and (entry['POS'] > self.stop + dist or entry['POS'] < self.start - dist): return False
        self.sig_snps.add(entry['SNP'])
        self.snps.update(set(entry['SP2'].replace('(1)','').split(',')))
        self.phenotypes.add((entry['group'], entry['pheno']))
        self.start = min(self.start, entry['POS'])
        self.stop = max(self.stop, entry['POS'])
        if entry['P'] < self.pval:
            self.pval = entry['P']
            self.sentinel_snp = entry['SNP']
        return True

    def to_dataframe(self):
        df = pd.DataFrame(dict(
            CHR = [self.chr], START = [self.start], STOP = [self.stop],
            SNP = [self.sentinel_snp], P = [self.pval],
            SIG = [','.join(sorted(self.sig_snps))]))
        for g, p in self.phenotypes: df[f'{g}_{p}'] = 1
        df['TOTAL'] = len(self.phenotypes)
        return df

def overlap_clumps(df, dist:int = 0):
    '''Set distance = 0 to merge clumps by SP2; otherwise merge clumps within specified distance'''
    df = df.sort_values(by = ['CHR','POS']).dropna()
    clumps = []
    current_clump = clump(df.iloc[0,:])
    for i in range(1, df.shape[0]):
        if not current_clump.update(df.iloc[i,:], dist):
            clumps.append(current_clump)
            current_clump = clump(df.iloc[i,:])
    clumps.append(current_clump)
    return pd.concat([c.to_dataframe() for c in clumps]).fillna(0), clumps

def parse_clump(pheno, clump_dir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/clump', pval = 5e-08):
    '''Identifies all clumps from a list of GWAS summary statistics'''
    # long format list of phenotypes
    if isinstance(pheno[0][1], list): pheno = [(g,p) for g, ps in pheno for p in ps]
    
    from _utils.path import find_clump
    clumps = []
    for g, p in pheno:
        try:
            clump_file,_ = find_clump(g, p, dirname = clump_dir, pval = pval)
            df = parse_clump_file(clump_file)
            df.insert(loc = 0, column = 'group', value = g)
            df.insert(loc = 1, column = 'pheno', value = p)
            clumps.append(df)
        except:
            Warning(f'No clump file found for {g}/{p}')
            continue
    clumps = pd.concat(clumps).sort_values(by = ['CHR','POS']).reset_index(drop = True)
    overlaps, _ = overlap_clumps(clumps)
    return clumps, overlaps
