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

import os
import pandas as pd
import numpy as np
from tqdm import tqdm
from fnmatch import fnmatch
from _utils.path import normaliser, find_gwas
from _utils.genetools import locus_to_name
from _utils.plugins.enrichr import enrichr_list, enrichr_to_revigo
from _utils.plugins.inrich import inrich
from multiprocessing import Pool, cpu_count

def parse_hyprcoloc_tabular(file):
    df = pd.read_table(file).drop(['iteration','dropped_trait'], axis=1) # valid clusters are now devoid of NA
    df.dropna(inplace = True)
    df = df.loc[df['regional_prob']>0.6, :] # filter by regional probability at 0.6 as in Nat Gen 2023
    return df

def parse_hyprcoloc_cluster(entry, pheno):
    '''
    entry: pandas DataFrame single row, with following columns: 'traits', 'candidate_snp';
    groups: a list output from _utils.path.find_gwas(long = False)
    '''
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

def _enrichr(df, traits):
    gene_list = locus_to_name(df, chrom_col = 'chromosome', start_col = 'start', stop_col = 'end')
    enrichr_result = enrichr_list(gene_list)
    return enrichr_result.assign(traits = traits)

def _inrich(df, traits):
    main, igt = inrich(df, chrom_col = 'chromosome', start_col = 'start', stop_col = 'end', niter = 10000)
    return main.assign(traits = traits), igt.assign(traits = traits)

def main(args):
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
    
    # enrichr and revigo analysis
    enrichr_res = []; revigo_res = [] 
    with Pool(processes = max(cpu_count()*2, orig.traits.unique().size)) as pool:
        enrichr_res = list(tqdm(pool.starmap(_enrichr, 
            [(group_df, phen_group) for phen_group, group_df in orig.groupby('traits')]), 
            total = orig.traits.unique().size))
        inrich_res = list(tqdm(pool.starmap(_inrich, 
            [(group_df, phen_group) for phen_group, group_df in orig.groupby('traits')]), 
            total = orig.traits.unique().size))
    revigo_res = enrichr_to_revigo(enrichr_res, keepcols = ['traits'])
    inrich_main = [x[0] for x in inrich_res]; inrich_igt = [x[1] for x in inrich_res]
    enrichr_res = pd.concat(enrichr_res)
    revigo_res = pd.concat(revigo_res)
    inrich_main = pd.concat(inrich_main)
    inrich_igt = pd.concat(inrich_igt)

    norm = normaliser()
    norm.normalise(orig).to_csv(f'{args.out}/{pheno_str}_coloc_raw.txt', sep = '\t', index = True)
    norm.normalise(summary).to_csv(f'{args.out}/{pheno_str}_coloc_summary.txt', sep = '\t', index = True)
    norm.normalise(clusters).to_csv(f'{args.out}/{pheno_str}_coloc_clusters.txt', sep = '\t', index = True)
    norm.normalise(enrichr_res).to_csv(f'{args.out}/{pheno_str}_coloc_enrichr.txt', sep = '\t', index = False)
    norm.normalise(revigo_res).to_csv(f'{args.out}/{pheno_str}_coloc_revigo.txt', sep = '\t', index = False)
    norm.normalise(inrich_main).to_csv(f'{args.out}/{pheno_str}_coloc_inrich_main.txt', sep = '\t', index = False)
    norm.normalise(inrich_igt).to_csv(f'{args.out}/{pheno_str}_coloc_inrich_igt.txt', sep = '\t', index = False)
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