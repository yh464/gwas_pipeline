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

# downstream analyses
def _stratify(df_score, adata, label):
    import pandas as pd
    df = adata.obs.loc[adata.obs.index.isin(df_score.index), :]
    cells = [df.index.tolist()]; annot = ['All']; cell_type = ['All']
    for l in label:
        if l not in df.columns: continue
        for group, df_group in df.groupby(l):
            cells.append(df_group.index.tolist())
            annot.append(l); cell_type.append(group)
    return pd.DataFrame(dict(cell = cells), index = pd.MultiIndex.from_frame(pd.DataFrame(dict(annot = annot, cell_type = cell_type))))

def _corr_pseudotime(df_score, adata, strata):
    import pandas as pd
    import numpy as np
    from tqdm import tqdm
    out_df = pd.DataFrame(columns = ['pseudotime'], index = strata.index, data = np.nan)

    for i, row in tqdm(strata.iterrows(), desc = 'Correlating scDRS score with pseudotime'):
        df_stratum = df_score.loc[row['cell'],:].join(adata.obs.loc[row['cell'], 'pseudotime'], how = 'inner').dropna()
        if df_stratum.shape[0] == 0: continue
        out_df.loc[i, 'pseudotime'] = df_stratum['norm_score'].corr(df_stratum['pseudotime'])
    return out_df

def _corr_genes(df_score, adata, strata, genes = [], block_size = 1000):
    import pandas as pd
    import numpy as np
    from tqdm import tqdm
    import gc
    if len(genes) == 0:
        if 'feature_type' in adata.var.columns:
            genes = adata.var_names[adata.var['feature_type'] == 'protein_coding'].tolist()
        elif 'highly_variable' in adata.var.columns:
            genes = adata.var_names[adata.var['highly_variable']].tolist()
        else: genes = adata.var_names.tolist()
    out_df = pd.DataFrame(columns = genes, index = strata.index, data = np.nan)

    # correlate in batches of 1000 genes (or manually specified block size)
    for gid in tqdm(range(0, len(genes), block_size), desc = f'Correlating scDRS with gene expression, total {int(len(genes)/block_size)+1} blocks'):
        gene_block = genes[gid:(gid+block_size)]
        expr = pd.DataFrame(data = adata[:, gene_block].to_memory().X.toarray(), index = adata.obs_names, columns = gene_block)
        for i, row in tqdm(strata.iterrows(), desc = f'Gene block {int(gid/block_size)+1}/{int(len(genes)/block_size)+1}'):
            df_stratum = df_score.loc[row['cell'],:].join(expr.loc[row['cell'], :], how = 'inner').dropna()
            if df_stratum.shape[0] == 0: continue
            out_df.loc[i, gene_block] = df_stratum[gene_block].corrwith(df_stratum['norm_score'])
        del expr; gc.collect()
    out_df = out_df.T
    med = out_df.median(axis = 1)
    out_df['med'] = med
    out_df = out_df.sort_values('med', ascending = False).drop(columns = 'med').T
    return out_df

def downstream_correlation(adata, df_score, label, genes = []):
    import pandas as pd
    strata = _stratify(df_score, adata, label)
    out_corr = []
    if 'pseudotime' in adata.obs.columns:
        out_corr.append(_corr_pseudotime(df_score, adata, strata))
    out_corr.append(_corr_genes(df_score, adata, strata, genes))
    return pd.concat(out_corr, axis = 1)

def main(args = None, **kwargs):
    from _utils.gadgets import namespace
    import os
    import scanpy as sc
    import scdrs
    import pandas as pd
    import matplotlib.pyplot as plt
    from _plots.aes import redblue_alpha
    from _plots import scatterplot_noaxis
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
        except: score = pd.read_table(out_score, index_col = 0, usecols = [0,2])
        if 'X_tsne' in adata.obsm.keys():
            x, y = adata.obsm['X_tsne'][:,0], adata.obsm['X_tsne'][:,1]
            rep = 'tSNE'
        elif 'X_umap' in adata.obsm.keys():
            x, y = adata.obsm['X_umap'][:,0], adata.obsm['X_umap'][:,1]
            rep = 'UMAP'
        else: raise ValueError('No tSNE or UMAP coordinates found in adata.obsm')
        scatterplot_noaxis(x, y, score['norm_score'], palette = redblue_alpha, s = 1, rep = rep)
        plt.savefig(out_fig, dpi = 400, bbox_inches = 'tight')
        plt.close()
        adata.file.close()
        if score.shape[1] < 3: del score

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
        adata.file.close()
    
    out_downstream = f'{args.out}.downstream.txt'
    if args.downstream and (not os.path.isfile(out_downstream) or args.force):
        adata = sc.read_h5ad(args.h5ad, 'r')
        score = pd.read_table(out_score, index_col = 0, usecols = [0, 2])
        corr = downstream_correlation(adata, score, args.label)
        corr.to_csv(out_downstream, index = True, sep = '\t')
        adata.file.close()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser('This script runs cell-type enrichments using scDRS')
    parser.add_argument('-i','--in', dest = '_in', help = 'Input gene list and weights', required = True)
    parser.add_argument('--h5ad', help = 'Input h5ad single-cell multiomics dataset', required = True)
    parser.add_argument('--label', nargs = '*', help = 'Columns containing cell classifications/types in the h5ad dataset',
        default = ['ROIGroup', 'ROIGroupCoarse', 'ROIGroupFine', 'roi', 'supercluster_term', 'cluster_id', 'subcluster_id', 'development_stage', # siletti
        'Class','Subclass','Type_updated', 'Cluster', 'Tissue']) # wang
    parser.add_argument('-d', '--downstream', help = 'Conduct downstream analyses', default = False, action = 'store_true')
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