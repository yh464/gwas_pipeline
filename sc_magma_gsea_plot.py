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

def main(args):
    import os
    import pandas as pd
    import scipy.stats as sts
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    from _plots import corr_heatmap
    from _utils.path import normaliser, find_gwas
    from tqdm import tqdm
    import warnings
    norm = normaliser()
    
    # find gene set files
    import os
    gsets = [(f'{args.gset}/{x[:-4]}', x[:-4]) for x in os.listdir(args.gset) if x[-4:] == '.txt']
    gscores = [(f'{args.gscore}/{x[:-4]}', x[:-4]) for x in os.listdir(args.gscore) if x[-4:] == '.txt']

    # identify phenotypes
    pheno = find_gwas(args.pheno, long = True)
    pheno_short = find_gwas(args.pheno)

    all_phenos = []
    for g, p in tqdm(pheno):
        all_gsets = []
        for gset_file, gset in gsets + gscores:
            magma_output = f'{args._in}/{g}/{p}/{p}.{args.annot}.{gset}.gsa.out'
            if not os.path.isfile(magma_output): warnings.warn(Warning(f'No MAGMA GSA output found for {g}/{p}/{gset}')); continue
            df = pd.read_table(magma_output, sep = '\\s+', comment = '#')
            df = df.rename(columns = {'FULL_NAME':'cell_type', 'P':'p', 'BETA':'beta'})
            if 'cell_type' not in df.columns: df['cell_type'] = df.VARIABLE
            df['cell_type'] = df['cell_type'].fillna(df.VARIABLE)

            if (gset_file,gset) in gsets: df.insert(0,'gene_set', value = gset)
            else: # gene score columns are named: <method>.<annotation>.<cell_type>
                gset_col = gset + '.' + df.cell_type.str.split('.', expand = True).iloc[:,1]
                df.insert(0,'gene_set', value = gset_col)

            # for 'clusters' and 'subclusters', map to the respective cell type term
            if os.path.isfile(f'{gset_file}.label'):
                labels = pd.read_table(f'{gset_file}.label', dtype = str)
                if 'label' in df.columns: df = df.rename(columns = {'label':'_orig_label'})
                df = df.merge(labels, how = 'left', on = 'cell_type')
            else: df['label'] = df['gene_set']
            df['label'] = df['label'].fillna(df['gene_set'])

            if not (gset_file,gset) in gsets:
                df.loc[:,'cell_type'] = ['.'.join(x.split('.')[2:]) for x in df.cell_type]
            df.insert(0,'pheno', value = p)
            df.insert(0,'group', value = g)

            # FDR correction
            q = df.p.values.copy()
            q = sts.false_discovery_control(q)
            df['q']= q
            all_gsets.append(df)
        if len(all_gsets) == 0: continue
        all_gsets = pd.concat(all_gsets, axis = 0)
        all_gsets.to_csv(f'{args._in}/{g}/{p}.{args.annot}.enrichments.txt', sep = '\t', index = False)
        all_phenos.append(all_gsets)

    if len(all_phenos) == 0: return
    all_phenos = pd.concat(all_phenos, axis = 0)
    out_prefix = f'{args._in}/'+'_'.join([x[0] for x in pheno_short]) + f'.{args.annot}.enrichments'
    all_phenos = norm.normalise(all_phenos)
    all_phenos.to_csv(f'{out_prefix}.txt', sep = '\t', index = False)

    # miami-like bar plot
    for gset in tqdm(all_phenos.gene_set.unique()):
      tmp = all_phenos.loc[all_phenos.gene_set == gset,:]
      if tmp.cell_type.unique().size < 1: continue
      tmp = tmp.assign(**{'-log(fdr)': (-np.log10(tmp.q) * (tmp['beta'] > 0))})
      if tmp.cell_type.unique().size < 100:
        fig1, ax1 = plt.subplots(tmp.cell_type.unique().size, 1, figsize = (len(pheno), 3*tmp.cell_type.unique().size), sharex = True, squeeze = False)
        ax1 = ax1[::-1,0] # invert y axis
        for i, ct in enumerate(tmp.cell_type.unique()):
            sns.barplot(tmp.loc[tmp.cell_type == ct, :], x = 'phenotype', y = '-log(fdr)', hue = 'group', ax = ax1[i], legend = False)
            if i > 0: ax1[i].set_xlabel('')
            ax1[i].axhline(-np.log10(0.05), color = 'k')
            ax1[i].axhline(np.log10(0.05), color = 'k')
        fig1.savefig(f'{out_prefix}.{gset}.barplot.pdf', bbox_inches = 'tight')
        plt.close(fig1)
      if tmp.cell_type.unique().size < 500: 
        fig = corr_heatmap(tmp[['group','phenotype','label','cell_type','beta','p','q']])
        fig.savefig(f'{out_prefix}.{gset}.pdf', bbox_inches = 'tight')
        plt.close(fig)
      
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'Parses MAGMA GSA outputs for a group of phenotypes')
    parser.add_argument('pheno', nargs = '*', help = 'Phenotypes')
    parser.add_argument('-i','--in', dest = '_in', help = 'MAGMA output directory',
      default = '../sc/magma_gsea')
    parser.add_argument('--annot', help = 'Annotation used to generate gene-level sumstats', default = 'ENSG')
    parser.add_argument('--gset', dest = 'gset', help = 'Gene sets to study enrichment, scans directory',
        default = '../multiomics/gene_set')
    parser.add_argument('--gscore', help = 'Directory containing gene scores', default = '../multiomics/gene_score')
    # no need for 'force'
    args = parser.parse_args()
    # path normalisation
    args.pheno.sort()
    import os
    for arg in ['_in','gset', 'gscore']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_var('%gset',r'.+','gene set')
    proj.add_input(args._in+'/%pheng/%pheno_%maf.%gset.gsa.out',__file__)
    proj.add_output(args._in+'/%pheng/%pheno_%maf.gsasummary.txt',__file__)
    try: main(args)
    except: cmdhistory.errlog()