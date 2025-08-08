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

def main(args):
    import os
    from _utils.path import normaliser, find_gwas
    import pandas as pd
    from fnmatch import fnmatch
    from logparser import crosscorr_parse
    
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
    else: 
        fout = f'{args.out}/corr_' + '_'.join(args.p1)
        summary_rev = summary.copy()
        summary_rev[['group1','pheno1','group2','pheno2']] = summary[['group2','pheno2','group1','pheno1']]
        summary = pd.concat([summary, summary_rev]).sort_values(['group1','pheno1','group2','pheno2']).drop_duplicates()
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