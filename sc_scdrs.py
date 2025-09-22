#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2025-08-21

Conducts scDRS scoring for a single phenotype

Requires following inputs: 
    Gene-level summary statistics
        Column specifications: gene ID (ENSGxxxxxxxx), p-value and effect size (optional)
    P-value threshold for sig genes (default fdr < 0.05)
    Minimum and maximum number of genes to include as sig genes (default: min 100, max 10% of genes)
    Pre-processed scDRS single cell data
        Column specifications: Columns of adata.obs that contain cell group labels
            (for Siletti et al 2023: ROIGroup, ROIGroupCoarse, ROIGroupFine, roi, supercluster_term, cluster_id, subcluster_id, development_stage)
'''

import matplotlib as mpl

def main(args = None, **kwargs):
    from _utils.gadgets import namespace
    import os
    import scanpy as sc
    import scdrs
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from _plots.aes import redblue
    from _plots.umap_scatterplot import scatterplot_noaxis
    import warnings
    if args == None:
        from _utils.gadgets import namespace
        args = namespace(**kwargs)

    # Stage 1: generate cell-specific scores
    out_score = f'{args.out}.score.txt'
    if not os.path.isfile(out_score) or args.force:
        weights = pd.read_table(args._in)
        gene_list = weights['gene'].tolist()
        gene_weight = weights['weight'].values if 'weight' in weights.columns else None
        adata = sc.read_h5ad(args.h5ad, 'r')
        block_size = int(2e9/adata.shape[1])
        if adata.shape[0] >= block_size:
            score_df = []
            for i in range(0, adata.shape[0], block_size):
                # control genes are the same as the same random seed is used
                cmin = i; cmax = min(i+block_size, adata.shape[0])
                tempfile = f'{args.out}.score.{cmin}_{cmax}.tmp.gz'
                if not os.path.isfile(tempfile) or args.force:
                    print(f'Scoring cells {cmin} - {cmax}')
                    temp = adata[cmin:cmax,:]
                    temp_score = scdrs.score_cell(temp, gene_list, gene_weight, return_ctrl_norm_score = True, verbose = True)
                    temp_score.to_csv(tempfile, index = True, sep = '\t')
                else: temp_score = pd.read_table(tempfile, index_col = 0)
                score_df.append(temp_score)
            score = pd.concat(score_df)
        else: score = scdrs.score_cell(adata, gene_list, gene_weight, return_ctrl_norm_score = True, verbose = True)
        score.to_csv(out_score, index = True, sep = '\t')
        adata.file.close()
    
    # Stage 2: plot cell_specific scores
    out_fig = f'{args.out}.score.png'
    if not os.path.isfile(out_fig) or args.force:
        adata = sc.read_h5ad(args.h5ad, 'r')
        try: score
        except: score = pd.read_table(out_score, index_col = 0)
        if 'X_tsne' in adata.obsm.keys():
            x, y = adata.obsm['X_tsne'][:,0], adata.obsm['X_tsne'][:,1]
            rep = 'tSNE'
        elif 'X_umap' in adata.obsm.keys():
            x, y = adata.obsm['X_umap'][:,0], adata.obsm['X_umap'][:,1]
            rep = 'UMAP'
        else: raise ValueError('No tSNE or UMAP coordinates found in adata.obsm')
        scatterplot_noaxis(x, y, score['norm_score'], palette = redblue, s = 0.1, rep = rep)
        # sns.set_theme(style = 'ticks
        # temp_df = pd.DataFrame(dict(tsne1 = adata.obsm['X_tsne'][:,0], tsne2 = adata.obsm['X_tsne'][:,1], score = score['norm_score']))
        # try: mpl.colormaps.register(redblue)
        # except: pass
        # sns.set_theme(style = 'ticks')
        # _, ax = plt.subplots(figsize = (5,5))
        # sns.scatterplot(temp_df, x = 'tsne1', y = 'tsne2', hue = 'score', palette = 'redblue', s = 1, ax = ax, edgecolor = None, linewidth = 0, legend = False)
        plt.savefig(out_fig, dpi = 400, bbox_inches = 'tight')
        plt.close()
        adata.file.close()

    # Stage 3: for each column of annotations, generate group-specific enrichments
    out_enrichment = f'{args.out}.enrichment.txt'
    if not os.path.isfile(out_enrichment) or args.force:
        adata = sc.read_h5ad(args.h5ad, 'r')
        try: score
        except: score = pd.read_table(out_score, index_col = 0)
        
        class_cols = [x for x in args.label if x in adata.obs.columns]
        mis_cols = [x for x in args.label if not x in adata.obs.columns]
        if len(mis_cols) > 0: warnings.warn('Following columns are missing from the h5ad dataset: '+' '.join(mis_cols))
        res = scdrs.method.downstream_group_analysis(adata, score, class_cols)
        enrichments = []
        for group, df in res.items():
            df.insert(loc = 0, column = 'cell_classification', value = group)
            enrichments.append(df)
        enrichments = pd.concat(enrichments)
        enrichments.to_csv(out_enrichment, index = True, sep = '\t')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser('This script runs cell-type enrichments using scDRS')
    parser.add_argument('-i','--in', dest = '_in', help = 'Input gene list and weights', required = True)
    parser.add_argument('--h5ad', help = 'Input h5ad single-cell multiomics dataset', required = True)
    parser.add_argument('--label', nargs = '*', help = 'Columns containing cell classifications/types in the h5ad dataset',
        default = ['ROIGroup', 'ROIGroupCoarse', 'ROIGroupFine', 'roi', 'supercluster_term', 'cluster_id', 'subcluster_id', 'development_stage', # siletti
        'Class','Subclass','Type_updated', 'Cluster', 'Tissue']) # wang
    parser.add_argument('-o', '--out', dest = 'out', help = 'output prefix, no .txt', required = True)
    parser.add_argument('-f','--force',dest = 'force', help = 'force overwrite', default = False, action = 'store_true')
    args = parser.parse_args()
    
    # path normalisation
    import os
    for arg in ['_in','h5ad','out']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))

    from _utils import cmdhistory
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()