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

def parse_hyprcoloc_tabular(file):
    import pandas as pd
    df = pd.read_table(file).drop(['iteration','dropped_trait'], axis=1) # valid clusters are now devoid of NA
    df.dropna(inplace = True)
    df = df.loc[df['regional_prob']>0.6, :] # filter by regional probability at 0.6 as in Nat Gen 2023
    return df

def parse_hyprcoloc_cluster(entry, pheno):
    '''
    entry: pandas DataFrame single row, with following columns: 'traits', 'candidate_snp';
    groups: a list output from _utils.path.find_gwas(long = False)
    '''
    import pandas as pd
    import numpy as np
    
    # initialise output
    snp = entry['candidate_snp']
    clusters = pd.DataFrame(data = np.nan, index = [snp], columns = [f'{g}_{p}' for g, ps in pheno for p in ps])
    summary = pd.DataFrame(data = np.nan, index = [snp], columns = [g for g,_ in pheno])
    
    for trait in entry['traits'].split(', '):
        clusters.loc[snp, trait.replace('/','_')] = 1
        summary.loc[snp, trait.split('/')[0]] = 1
    clusters.insert(loc = 0, column = 'regional_prob', value = entry['regional_prob'])
    clusters.insert(loc = 1, column = 'prob_explained_by_snp', value = entry['posterior_explained_by_snp'])
    summary.insert(loc = 0, column = 'regional_prob', value = entry['regional_prob'])
    summary.insert(loc = 1, column = 'prob_explained_by_snp', value = entry['posterior_explained_by_snp'])
    return clusters, summary
    
    
def main(args):
    import os
    import pandas as pd
    from fnmatch import fnmatch
    from _utils.path import normaliser, find_gwas
    
    pheno = find_gwas(args.pheno, dirname = args.gwa, long = False)
    pheno_str = '_'.join([g for g,_ in pheno])
    in_dir = f'{args._in}/{pheno_str}'  
    
    orig = []
    summary = []
    clusters = []
    for x in sorted(os.listdir(in_dir)):
        if not fnmatch(x, '*hyprcoloc.txt'): continue
        # parse chromosomal position
        tmp = x.split('_')
        chrom = int(tmp[0][3:]); start = int(tmp[1]); end = int(tmp[2])
        
        df = parse_hyprcoloc_tabular(f'{in_dir}/{x}')
        for i in range(df.shape[0]):
            c, s = parse_hyprcoloc_cluster(df.iloc[i,:], pheno)
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
    
    norm = normaliser()
    norm.normalise(orig).to_csv(f'{args.out}/{pheno_str}_coloc_raw.txt', sep = '\t', index = True)
    norm.normalise(summary).to_csv(f'{args.out}/{pheno_str}_coloc_summary.txt', sep = '\t', index = True)
    norm.normalise(clusters).to_csv(f'{args.out}/{pheno_str}_coloc_clusters.txt', sep = '\t', index = True)
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