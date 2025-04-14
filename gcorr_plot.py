#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2023-07-18
Version 2: 2024-11-14
Version 3: 2025-04-09

A simplified script to plot genetic correlation between groups of phenotypes

Requires following inputs: 
    LDSC rg logs
Changelog:
    Changed the file name structure (now scans sumstats directory to find relevant rg log files)
    Changed the plotting style to scatterplot-style heatmaps
    Applied a wide- and long-format tabular output
'''

def crosscorr_parse(gwa1, gwa2 = [], 
        logdir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/gcorr/rglog',
        h2dir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/gcorr/ldsc_sumstats',
        exclude = []):
    '''
    gwa1 and gwa2 are lists of (group, pheno_list) tuples from logparser.find_gwas(long=False)
    leave gwa2 blank to estimate auto-correlations of gwa1
    '''
    import os
    import numpy as np
    import pandas as pd
    import scipy.stats as sts
    from logparser import parse_rg_log, parse_h2_log
    summary = []
    
    from _utils.path import pair_gwas
    pairwise = pair_gwas(gwa1, gwa2)
    
    for g1, p1s, g2, p2s in pairwise:
        if g1 > g2: 
            g1, p1s, g2, p2s = g2, p2s, g1, p1s
            flip = True
        for p1 in p1s:
            fname = f'{logdir}/{g1}.{g2}/{g1}_{p1}.{g2}.rg.log'
            if not os.path.isfile(fname): continue
            rg = parse_rg_log(fname)
            rg['fixed_int'] = False
            if flip: rg.iloc[:,[0,1,2,3]] = rg.iloc[:,[2,3,0,1]]
            summary.append(rg)
            
            fname_noint = fname.replace('.rg.log','.noint.rg.log')
            if os.path.isfile(fname_noint):
                rg = parse_rg_log(fname_noint)
                rg['fixed_int'] = True
                if flip: rg.iloc[:,[0,1,2,3]] = rg.iloc[:,[2,3,0,1]]
                summary.append(rg)
            
            if g1 == g2 and h2dir != None: # heritability
                fname = f'{h2dir}/{g1}/{p1}.h2.log'
                h2, se = parse_h2_log(fname)
                rg = pd.DataFrame(dict(group1 = [g1], pheno1 = p1, group2 = g1, 
                    pheno2 = p1, rg = h2, se = se, 
                    p = 1-sts.chi2.cdf(rg**2/se**2, df = 1),fixed_int = False))
                if flip: rg.iloc[:,[0,1,2,3]] = rg.iloc[:,[2,3,0,1]]
                summary.append(rg)
            
    summary = pd.concat(summary) # creates a long format table
    summary.insert(loc = len(summary.columns), column = 'q', value = np.nan)
    for g1,p1s in gwa1: # FDR correction for each IDP, which are non-independent
        for p1 in p1s:
            summary.loc[(summary.pheno1==p1) & (summary.group1==g1) & ~np.isnan(summary.p),'q'] \
                = sts.false_discovery_control(
            summary.loc[(summary.pheno1==p1)&(summary.group1==g1) & ~np.isnan(summary.p),'p'])
    
    return summary

def main(args):
    import os
    from _utils.path import normaliser, find_gwas
    from fnmatch import fnmatch
    
    # scans directories to include sumstats
    gwa1 = find_gwas(*args.p1, dirname = args.sumstats, ext = 'sumstats')
    gwa2 = find_gwas(*args.p2, dirname = args.sumstats, ext = 'sumstats')

    # parse LDSC log files
    summary = crosscorr_parse(gwa1, gwa2, args._in, h2dir = args.sumstats, exclude = args.exclude)
    
    # drop --exclude phenotypes in summary
    to_drop = []
    for idx in summary.index:
        if any([fnmatch(summary.loc[idx,'pheno1'], x) for x in args.exclude]) or \
            any([fnmatch(summary.loc[idx,'pheno2'], x) for x in args.exclude]):
            to_drop.append(idx)
    summary = summary.drop(to_drop, axis = 0)
    
    # tabular output, wide and long
    norm = normaliser()
    if len(args.p2) > 0:
        fout = f'{args.out}/crosscorr_' + '_'.join(args.p1)+'.'+'_'.join(args.p2)
    else: fout = f'{args.out}/corr_' + '_'.join(args.p1)
    rg_tbl = summary.pivot_table(index = ['group1','pheno1'], columns = ['group2','pheno2'], values = 'rg')
    norm.normalise(rg_tbl).to_csv(f'{fout}.wide.txt', index_label = False, sep = '\t',
                  header = True, index = True)
    summary = norm.normalise(summary)
    if not os.path.isfile(f'{fout}.txt') or len(args.exclude) == 0:
        summary.to_csv(f'{fout}.txt', index = False, sep = '\t')
    
    # plot figure
    from _plots import corr_heatmap
    fig = corr_heatmap(summary, annot = 'Heritability')
    fig.savefig(f'{fout}.pdf', bbox_inches = 'tight')
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'This programme parses genetic cross-correlation between two groups of phenotypes')
    parser.add_argument('-p1', help = 'First group of phenotypes to correlate', nargs = '*',
      default = [])
    parser.add_argument('-p2', help = 'Second group of phenotypes to correlate', nargs = '*',
      default = [])
    parser.add_argument('-i','--in', dest = '_in', help = 'Directory for rg logs',
      default = '../gcorr/rglog/')
    parser.add_argument('--exclude', help = 'phenotypes to exclude', nargs = '*', default = [])
    parser.add_argument('--sumstats', help = 'sumstats directory to be scanned for file names',
      default = '../gcorr/ldsc_sumstats/')
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory')
    # always forces output
    args = parser.parse_args()
    
    import os
    args.p1.sort()
    args.p2.sort()
    args._in = os.path.realpath(args._in)
    args.sumstats = os.path.realpath(args.sumstats)
    if type(args.out) == type(None): args.out = os.path.realpath(f'{args._in}/../')
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheno.%pheno.rg.log', __file__)
    proj.add_output(args.out+'/crosscorr_.*.pdf', __file__)
    try: main(args)
    except: cmdhistory.errlog()