#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2023-07-20
Version 2: 2024-11-14
Version 3: 2025-09-11

Summarises gene-set level enrichment for HMAGMA and MAGMA outputs

Preceding workflow:
    annot_batch.py
Requires following inputs:
    MAGMA GSA outputs
Changelog:
    Added the 'annot' paramter to reflect the multitude of Hi-C and nearest-gene-based annotations
    Changed the heatmap to a scatterplot-style heatmap
'''

from ast import parse
import os
from fnmatch import fnmatch
import pandas as pd
import scipy.stats as sts
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from _plots import corr_heatmap
from _utils.path import normaliser, find_gwas
from tqdm import tqdm
from multiprocessing import Pool
import warnings
from _utils import logger
log = logger.logger()

def parse_gset_name(df, gset, gset_file):
    gset_df = df.cell_type.str.split('.', expand = True)
    if gset_df.shape[1] == 1: 
        df.insert(0,'gene_set', value = gset) # no annotation
        df['cell_type_parsed'] = df['cell_type']
    elif gset_df.shape[0] == 2 and gset_df.iloc[:,0].unique().size < gset_df.iloc[:,1].unique().size: 
        df.insert(0,'gene_set', value = gset + '.' + gset_df.iloc[:,0]) # gene sets are annotated as <annotation>.<cell_type>
        df['cell_type_parsed'] = gset_df.iloc[:,1] # cell type is in the second column
    elif gset_df.shape[1] == 3 and gset_df.iloc[:,0].unique().size == 1:
        df.insert(0,'gene_set', value = gset + '.' + gset_df.iloc[:,1]) # gene sets are annotated as <method>.<annotation>.<cell_type>
        df['cell_type_parsed'] = gset_df.iloc[:,2] # cell type is in the third column
    else: 
        df.insert(0,'gene_set', value = gset) # fallback to no annotation
        df['cell_type_parsed'] = df['cell_type']

    if os.path.isfile(f'{gset_file}.label'):
        labels = pd.read_table(f'{gset_file}.label', dtype = str)
        if 'label' in df.columns: df = df.rename(columns = {'label':'_orig_label'})
        df = df.merge(labels, how = 'left', on = 'cell_type')
    else: df['label'] = df['gene_set']
    df['label'] = df['label'].fillna(df['gene_set'])
    df = df.drop(columns = ['cell_type']).rename(columns = {'cell_type_parsed':'cell_type'})
    return df

def process_pheno(gsets, g, p, args):
    all_gsets = []
    most_sig = []
    for gset_file, gset in gsets:
        # process MAGMA output for main analysis
        magma_output = f'{args._in}/{g}/{p}/{p}.{args.annot}.{gset}.gsa.out'
        if not os.path.isfile(magma_output): log.warn(Warning(f'No MAGMA GSA output found for {g}/{p}/{gset}')); continue
        df = pd.read_table(magma_output, sep = '\\s+', comment = '#')
        df = df.rename(columns = {'FULL_NAME':'cell_type', 'P':'p', 'BETA_STD':'beta', 'BETA': 'beta_raw'})
        if 'cell_type' not in df.columns: df['cell_type'] = df.VARIABLE
        df['cell_type'] = df['cell_type'].fillna(df.VARIABLE)

        # parse gene set names
        df = parse_gset_name(df, gset, gset_file)
        df.insert(0,'phenotype', value = p)
        df.insert(0,'group', value = g)

        # FDR correction
        q = df.p.values.copy()
        q = sts.false_discovery_control(q)
        df['q']= q
        all_gsets.append(df)
        df = df.sort_values('p').reset_index(drop = True)
        most_sig.append(df.loc[df.q < 0.05,:])

        # process output of conditional analysis
        cond_output = f'{args._in}/{g}/{p}/{p}.{args.annot}.{gset}.cond.gsa.out'
        if not os.path.isfile(cond_output): continue
        df = pd.read_table(cond_output, sep = '\\s+', comment = '#')
        df = df.rename(columns = {'VARIABLE':'cell_type', 'P':'p', 'BETA_STD':'beta', 'BETA': 'beta_raw'})
        df = parse_gset_name(df, gset, gset_file)
        df['analysed_cell_type'] = df['cell_type']
        # reorder cell types 1, 0, 3, 2, 5, 4, ... to get the cell types being conditioned on
        reorder = np.arange(df.shape[0])
        reorder = np.stack([reorder[1::2], reorder[0::2]]).T.reshape(-1) 
        df['conditioned_on'] = df['cell_type'].iloc[reorder].values
        df['tmp1'] = 'conditioned_on'; df['tmp2'] = 'analysed_cell_type'
        # FDR correction with respect to each trait being conditioned on
        fig = corr_heatmap(df[['tmp1','conditioned_on','tmp2','analysed_cell_type','beta','p']])
        fig.savefig(f'{args._in}/{g}/{p}.{args.annot}.{gset}.cond.heatmap.pdf', bbox_inches = 'tight')
        plt.close(fig)

    most_sig = pd.concat(most_sig, axis = 0) if len(most_sig) > 0 else None
    if len(all_gsets) == 0: return None
    all_gsets = pd.concat(all_gsets, axis = 0)
    normaliser().normalise(all_gsets).to_csv(f'{args._in}/{g}/{p}.{args.annot}.enrichments.txt', sep = '\t', index = False)
    return all_gsets, most_sig

@log.profile
def main(args):   
    # find gene set files
    gsets = [(f'{args.gset}/{x[:-4]}', x[:-4]) for x in os.listdir(args.gset) if x[-4:] == '.txt']
    gscores = [(f'{args.gscore}/{x[:-4]}', x[:-4]) for x in os.listdir(args.gscore) if x[-4:] == '.txt']
    if len(args.subset) > 0:
        gsets = [x for x in gsets if any([fnmatch(x[1], f'*{s}*') for s in args.subset])]
        gscores = [x for x in gscores if any([fnmatch(x[1], f'*{s}*') for s in args.subset])]
    log.log(f'Found {len(gsets)} gene sets and {len(gscores)} gene scores to process')
    for _, x in gsets + gscores:
        log.log(f'    {x}')

    # identify phenotypes
    pheno = find_gwas(args.pheno, long = True)
    pheno_short = find_gwas(args.pheno)

    with Pool(processes = min(8, len(pheno))) as pool:
        results = list(tqdm(pool.starmap(process_pheno, [(gsets + gscores, g, p, args) for g, p in pheno]), total = len(pheno)))
    all_phenos = [x[0] for x in results if x is not None]
    most_sig = [x[1] for x in results if x is not None]

    if len(most_sig) == 0: return
    most_sig = pd.concat(most_sig, axis = 0)
    most_sig.to_clipboard(sep = '\t', index = False)
    most_sig = most_sig.sort_values('p').reset_index(drop = True)
    log.log(most_sig.head())

    if len(all_phenos) == 0: return
    all_phenos = pd.concat(all_phenos, axis = 0)
    out_prefix = f'{args._in}/{args.annot}'
    # all_phenos.to_csv(f'{out_prefix}.txt', sep = '\t', index = False)

    # miami-like bar plot
    for gset, gset_df in tqdm(all_phenos.groupby('gene_set'), total = all_phenos.gene_set.unique().size):
      if gset.find('GO') >= 0: # for GO terms, just report significant terms in a table
        gset_df = gset_df.loc[gset_df.p < 0.05,:]
        gset_df = gset_df.sort_values('p').reset_index(drop = True)
        gset_df.to_csv(f'{out_prefix}.{gset}.txt', sep = '\t', index = False)
        continue

      gset_df.to_csv(f'{out_prefix}.{gset}.txt', sep = '\t', index = False)
      if gset_df.cell_type.unique().size < 1: continue
      gset_df = gset_df.assign(**{'-log(fdr)': (-np.log10(gset_df.q) * (gset_df['beta'] > 0))})
      if gset_df.cell_type.unique().size < 100:
        fig1, ax1 = plt.subplots(gset_df.cell_type.unique().size, 1, figsize = (len(pheno), 3*gset_df.cell_type.unique().size), sharex = True, squeeze = False)
        ax1 = ax1[::-1,0] # invert y axis
        for i, ct in enumerate(gset_df.cell_type.unique()):
            sns.barplot(gset_df.loc[gset_df.cell_type == ct, :], x = 'phenotype', y = '-log(fdr)', hue = 'group', ax = ax1[i], legend = False)
            if i > 0: ax1[i].set_xlabel('')
            ax1[i].axhline(-np.log10(0.05), color = 'k')
            ax1[i].axhline(np.log10(0.05), color = 'k')
        fig1.savefig(f'{out_prefix}.{gset}.barplot.pdf', bbox_inches = 'tight')
        plt.close(fig1)
      if gset_df.cell_type.unique().size < 500: 
        fig = corr_heatmap(gset_df[['group','phenotype','label','cell_type','beta','p','q']])
        fig.savefig(f'{out_prefix}.{gset}.pdf', bbox_inches = 'tight')
        plt.close(fig)
      
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'Parses MAGMA GSA outputs for a group of phenotypes')
    parser.add_argument('pheno', nargs = '*', help = 'Phenotypes')
    parser.add_argument('-i','--in', dest = '_in', help = 'MAGMA output directory',
      default = '../sc/magma_gsea')
    parser.add_argument('-s', '--subset', help = 'Subset of gene sets and gene scores to analyse', nargs = '*', default = [])
    parser.add_argument('-a', '--annot', help = 'Annotation used to generate gene-level sumstats', default = 'ENSG_10kb')
    parser.add_argument('--gset', dest = 'gset', help = 'Gene sets to study enrichment, scans directory',
        default = '../multiomics/gene_set')
    parser.add_argument('--gscore', help = 'Directory containing gene scores', default = '../multiomics/gene_score')
    # always overwrites
    args = parser.parse_args()
    # path normalisation
    args.pheno.sort()
    import os
    for arg in ['_in','gset', 'gscore']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    
    from _utils import cmdhistory, path
    log.splash(args)
    cmdhistory.log()
    proj = path.project()
    try: main(args)
    except: cmdhistory.errlog()