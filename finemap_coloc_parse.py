#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-01-08

A script to summarise hyprcoloc output tables

Preceding workflow:
    finemap_coloc_batch.py
Required input:
    hyprcoloc tabular output
'''

def parse_coloc_tabular(file):
    import pandas as pd
    df = pd.read_table(file).drop(['iteration','dropped_trait'], axis=1) # valid clusters are now devoid of NA
    df.dropna(inplace = True)
    df = df.loc[df['regional_prob']>0.6, :] # filter by regional probability at 0.6 as in Nat Gen 2023
    return df

def parse_coloc_cluster(entry, groups):
    '''
    entry: pandas DataFrame single row, with following columns: 'traits', 'candidate_snp';
    groups: pandas DataFrame with 'group' and 'pheno' columns
    '''
    import pandas as pd
    import numpy as np
    groups['names'] = groups['group'] + '_' + groups['pheno']
    
    # initialise output
    snp = entry['candidate_snp']
    clusters = pd.DataFrame(data = np.nan, index = [snp], columns = groups['names'])
    summary = pd.DataFrame(data = np.nan, index = [snp], columns = groups['group'].unique())
    
    for trait in entry['traits'].split(', '):
        group = groups.loc[groups.names==trait, 'group']
        clusters.loc[snp, trait] = 1
        summary.loc[snp, group] = 1
    clusters.insert(loc = 0, column = 'regional_prob', value = entry['regional_prob'])
    clusters.insert(loc = 1, column = 'prob_explained_by_snp', value = entry['posterior_explained_by_snp'])
    summary.insert(loc = 0, column = 'regional_prob', value = entry['regional_prob'])
    summary.insert(loc = 1, column = 'prob_explained_by_snp', value = entry['posterior_explained_by_snp'])
    return clusters, summary
    
    
def main(args):
    import os
    import pandas as pd
    from fnmatch import fnmatch
    from _utils.path import normaliser
    
    pheno = '_'.join(sorted(args.pheno))
    os.chdir(f'{args._in}/{pheno}')
    
    groups = []
    for p in args.pheno:
        for x in os.listdir(f'{args.gwa}/{p}'):
            if not fnmatch(x, '*.fastGWA'): continue
            if fnmatch(x, '*_X.fastGWA'): continue
            if fnmatch(x, '*_all_chrs*'): continue
            prefix = x.replace('.fastGWA','')
            groups.append(pd.DataFrame(dict(group = p, pheno = [prefix])))
    groups = pd.concat(groups).reset_index(drop = True)        
    
    orig = []
    summary = []
    clusters = []
    for x in sorted(os.listdir(f'{args._in}/{pheno}')):
        if not fnmatch(x, '*coloc.txt'): continue
        # parse chromosomal position
        tmp = x.split('_')
        chrom = int(tmp[0][3:]); start = int(tmp[1]); end = int(tmp[2])
        
        df = parse_coloc_tabular(x)
        for i in range(df.shape[0]):
            c, s = parse_coloc_cluster(df.iloc[i,:], groups)
            c.insert(loc = 0, column = 'CHR', value = chrom)
            c.insert(loc = 1, column = 'START', value = start)
            c.insert(loc = 2, column = 'END', value = end)
            
            s.insert(loc = 0, column = 'CHR', value = chrom)
            s.insert(loc = 1, column = 'START', value = start)
            s.insert(loc = 2, column = 'END', value = end)
            summary.append(s)
            clusters.append(c)
        df = df[['candidate_snp','regional_prob','posterior_prob','posterior_explained_by_snp','traits']]
        df = df.set_index('candidate_snp')
        df.insert(loc = 0, column = 'chromosome', value = chrom)
        df.insert(loc = 1, column = 'start', value = start)
        df.insert(loc = 2, column = 'end', value = end)
        orig.append(df)
    orig = pd.concat(orig)
    summary = pd.concat(summary)
    summary.sort_values(by = summary.columns.tolist(), inplace = True)
    clusters = pd.concat(clusters)
    summary.index.name = 'SNP'
    clusters.index.name = 'SNP'
    clusters = clusters.T
    
    norm = normaliser()
    norm.normalise(orig).to_csv(f'{args.out}/{pheno}_coloc_raw.txt', sep = '\t', index = True)
    norm.normalise(summary).to_csv(f'{args.out}/{pheno}_coloc_summary.txt', sep = '\t', index = True)
    norm.normalise(clusters).to_csv(f'{args.out}/{pheno}_coloc_clusters.txt', sep = '\t', index = True)
    return

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
      description = 'This programme summarises hyprcoloc output')
    parser.add_argument('pheno', help = 'Phenotype groups', nargs = '*')
    parser.add_argument('-i','--in', dest = '_in', help = 'Directory containing coloc outputs',
      default = '../coloc/')
    parser.add_argument('-g','--gwa', dest = 'gwa', help = 'Directory containing all summary stats (for phenotype names only)',
      default = '../gwa/')
    parser.add_argument('-o', '--out', dest = 'out', help = 'output directory')
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
      default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    args._in = os.path.realpath(args._in)
    if type(args.out) == type(None): args.out = args._in
    else: args.out = os.path.realpath(args.out)
    args.gwa = os.path.realpath(args.gwa)
    
    from _utils import path, cmdhistory, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_var('%pval',r'[0-9.+-e]+', 'minor allele freq') # only allows digits and decimals
    proj.add_input(args._in+'/%pheng/*coloc.txt', __file__)
    proj.add_output(args.out+'/%pheng_coloc_summary.txt',__file__)
    proj.add_output(args.out+'/%pheng_coloc_clusters.txt',__file__)
    try: main(args)
    except: cmdhistory.errlog()