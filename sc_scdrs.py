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

import pandas as pd
import numpy as np
from tqdm import tqdm
import scanpy as sc
import scdrs, gget, gc, warnings, os, argparse
from ._utils.genetools import ensg_to_name
from multiprocessing import Pool, cpu_count
from ._utils.gadgets import namespace
import matplotlib.pyplot as plt
from ._plots.aes import redblue_alpha
from ._plots import scatterplot_noaxis, temporal_regplot
from ._utils.gadgets import mv_symlink

# downstream analyses
def _stratify(df_score, adata, label):
    df = adata.obs.loc[adata.obs.index.isin(df_score.index), :]
    cells = [df.index.tolist()]; annot = ['All']; cell_type = ['All']
    for l in label:
        if l not in df.columns: continue
        for group, df_group in df.groupby(l):
            cells.append(df_group.index.tolist())
            annot.append(l); cell_type.append(group)
    return pd.DataFrame(dict(cell = cells), index = pd.MultiIndex.from_frame(pd.DataFrame(dict(annot = annot, cell_type = cell_type))))

def _corr_pseudotime(df_score, adata, strata):
    out_df = pd.DataFrame(columns = ['pseudotime'], index = strata.index, data = np.nan)

    for i, row in tqdm(strata.iterrows(), desc = 'Correlating scDRS score with pseudotime'):
        df_stratum = df_score.loc[row['cell'],:].join(adata.obs.loc[row['cell'], 'pseudotime'], how = 'inner').dropna()
        if df_stratum.shape[0] == 0: continue
        out_df.loc[i, 'pseudotime'] = df_stratum['norm_score'].corr(df_stratum['pseudotime'])
    return out_df

def _corr_genes(df_score, adata, strata, genes = [], block_size = 1000):
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
    strata = _stratify(df_score, adata, label)
    out_corr = []
    if 'pseudotime' in adata.obs.columns:
        out_corr.append(_corr_pseudotime(df_score, adata, strata))
    out_corr.append(_corr_genes(df_score, adata, strata, genes))
    return pd.concat(out_corr, axis = 1)

def _enrichr(stratum, databases = ['GO_Biological_Process_2025','SynGO_2024'], top = [100, 200, 500, 1000]):
    stratum = stratum.dropna().sort_values(ascending = True)
    gene_list = stratum.index.tolist()
    out = []
    for db in databases:
        for n in top:
            enrichr_res = gget.enrichr(gene_list[:n], db)
            out.append(pd.DataFrame(dict(
                annot = stratum.name[0], cell_type = stratum.name[1], n_genes = len(stratum), top = n, database = db,
                sign = '-', process = enrichr_res.loc[:20,'path_name'], p = enrichr_res.loc[:20,'p_val'])))
            enrichr_res = gget.enrichr(gene_list[-n:], db)
            out.append(pd.DataFrame(dict(
                annot = stratum.name[0], cell_type = stratum.name[1], n_genes = len(stratum), top = n, database = db,
                sign = '+', process = enrichr_res.loc[:20,'path_name'], p = enrichr_res.loc[:20,'p_val'])))
    return pd.concat(out, axis = 0)

def downstream_enrichr(corr_df):
    # only include classes with <=100 cell types
    strata = corr_df.index.to_frame()
    for a in strata['annot'].unique():
        if strata.loc[strata['annot'] == a,:].shape[0] > 100:
            strata = strata.loc[strata['annot'] != a,:]
    strata = pd.MultiIndex.from_frame(strata)
    corr_df = corr_df.loc[strata, corr_df.columns != 'pseudotime']
    corr_df.columns = ensg_to_name(corr_df.columns.tolist())

    with Pool(min(cpu_count()*4, len(strata))) as p:
        summary = list(tqdm(p.imap(_enrichr, [corr_df.loc[i,:] for i in strata]), total = len(strata), desc = 'Conducting Enrichr analysis'))
    enrichr_summary = pd.concat(summary, axis = 0)
    
    from ._plugins.enrichr import enrichr_to_revigo
    revigo_summary = enrichr_to_revigo(
        [df for _, df in enrichr_summary.groupby(['annot','cell_type','n_genes','top','sign'])],
        name_col = 'process', pval_col = 'p'
    )
    for idx, (group, _) in enumerate(enrichr_summary.groupby(['annot','cell_type','n_genes','top','sign'])):
        revigo_summary[idx] = revigo_summary[idx].assign(
            annot = group[0], cell_type = group[1], n_genes = group[2], top = group[3], sign = group[4]
        )
    revigo_summary = pd.concat(revigo_summary, axis = 0)
    return enrichr_summary, revigo_summary

def main(args = None, **kwargs):
    if args == None: args = namespace(**kwargs)

    # Stage 1: generate cell-specific scores
    out_score = f'{args.out}.score.txt'
    if not os.path.isfile(out_score) or args.force:
        weights = pd.read_table(args._in)
        gene_list = weights['gene'].iloc[:args.nsig].tolist()
        gene_weight = weights['weight'].iloc[:args.nsig].values if 'weight' in weights.columns else None
        adata = sc.read_h5ad(args.h5ad, 'r')
        block_size = int(2e9/adata.shape[1])
        if adata.shape[0] >= block_size:
            score_df = []
            for i in range(0, adata.shape[0], block_size):
                # control genes are the same as the same random seed is used
                cmin = i; cmax = min(i+block_size, adata.shape[0])
                tempfile = f'{args.out}.score.{cmin}_{cmax}.tmp.gz'
                try: 
                    if args.force: raise FileNotFoundError
                    temp_score = pd.read_table(tempfile, index_col = 0)
                except:
                    print(f'Scoring cells {cmin} - {cmax}')
                    temp = adata[cmin:cmax,:]
                    temp_score = scdrs.score_cell(temp, gene_list, gene_weight, return_ctrl_norm_score = True, verbose = True)
                    temp_score.to_csv(tempfile, index = True, sep = '\t')
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
    
    # Downstream analysis 1: for each cell type, correlate scDRS score with pseudotime and gene expression
    score = pd.read_table(out_score, index_col = 0, usecols = [0,2])
    adata = sc.read_h5ad(args.h5ad, 'r')
    out_downstream = f'{args.out}.downstream.txt'
    if args.downstream and (not os.path.isfile(out_downstream) or args.force):
        corr = downstream_correlation(adata, score, args.label)
        corr.to_csv(out_downstream, index = True, sep = '\t')
        adata.file.close()

    # Downstream analysis 2: for each cell type, conduct enrichment analysis using top correlated genes
    out_enrichr = f'{args.out}.downstream.enrichr.txt'
    out_revigo = f'{args.out}.downstream.revigo.txt'
    if args.downstream and (not os.path.isfile(out_enrichr) or not os.path.isfile(out_revigo) or args.force):
        try: corr
        except: corr = pd.read_table(out_downstream, index_col = [0,1])
        enrichr, revigo = downstream_enrichr(corr)
        enrichr.to_csv(out_enrichr, index = False, sep = '\t')
        revigo.to_csv(out_revigo, index = False, sep = '\t')

    # Downstream analysis 3: plot scDRS score with pseudotime, stratified by cell type
    out_pseudotime_fig = f'{args.out}.pseudotime.png'
    cell_type_cols = [x for x in args.label if x.lower().find('type') > -1 or x.lower().find('annot') > -1 and x in adata.obs.columns]
    if args.downstream and 'pseudotime' in adata.obs.columns and len(cell_type_cols) > 0 and \
        (not os.path.isfile(out_pseudotime_fig) or args.force):
        df = pd.concat([score['norm_score'], adata.obs[['pseudotime', cell_type_cols[0]]]], axis = 1).dropna()
        fig = temporal_regplot(
            df, x = 'pseudotime', y = 'norm_score', hue = cell_type_cols[0],
            xlabel = '', ylabel = ' scDRS score', order = 2, clip_tail = 0.025
        )
        fig.savefig(out_pseudotime_fig, dpi = 400, bbox_inches = 'tight')
        plt.close()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser('This script runs cell-type enrichments using scDRS')
    parser.add_argument('-i','--in', dest = '_in', help = 'Input gene list and weights', required = True)
    parser.add_argument('--h5ad', help = 'Input h5ad single-cell multiomics dataset', required = True)
    parser.add_argument('--label', nargs = '*', help = 'Columns containing cell classifications/types in the h5ad dataset',
        default = ['ROIGroup', 'ROIGroupCoarse', 'ROIGroupFine', 'roi', 'supercluster_term', 'cluster_id', 'subcluster_id', 'development_stage', # siletti
        'Class','Subclass','Type_updated', 'Cluster', 'Tissue']) # wang
    parser.add_argument('-n','--nsig', help = 'Number of significant genes, default 1000', 
        type = int, default = 1000)
    parser.add_argument('-d', '--downstream', help = 'Conduct downstream analyses', default = False, action = 'store_true')
    parser.add_argument('-o', '--out', dest = 'out', help = 'output prefix, no .txt', required = True)
    parser.add_argument('-f','--force',dest = 'force', help = 'force overwrite', default = False, action = 'store_true')
    args = parser.parse_args()
    
    # path normalisation
    import os
    for arg in ['_in','h5ad','out']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))

    from ._utils import cmdhistory
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()