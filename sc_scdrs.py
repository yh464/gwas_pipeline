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

def main(args = None, **kwargs):
    from _utils.gadgets import namespace
    import os
    import scdrs
    import pandas as pd
    if args == None:
        from _utils.gadgets import namespace
        args = namespace(**kwargs)

    # Stage 1: generate cell-specific scores
    out_score = f'{args.out}.score.txt'
    if not os.path.isfile(out_score) or args.force:
        weights = pd.read_table(args._in)
        gene_list = weights['gene'].tolist()
        gene_weight = weights['weight'].values if 'weight' in weights.columns else None
        adata = scdrs.util.load_h5ad(args.h5ad)
        score = scdrs.score_cell(adata, gene_list, gene_weight, return_ctrl_norm_score = True, verbose = True)
        score.to_csv(out_score, index = True, sep = '\t')
    
    # Stage 2: for each column of annotations, generate group-specific enrichments
    out_enrichment = f'{args.out}.enrichment.txt'
    if not os.path.isfile(out_enrichment) or args.force:
        try: adata
        except: adata = scdrs.util.load_h5ad(args.h5ad)
        try: score
        except: score = pd.read_table(out_score, index_col = 0)
        
        class_cols = [x for x in args.label if x in adata.obs.columns]
        mis_cols = [x for x in args.label if not x in adata.obs.columns]
        if len(mis_cols) > 0: Warning('Following columns are missing from the h5ad dataset: '+' '.join(mis_cols))
        res = scdrs.method.downstream_group_analysis(adata, score, class_cols)
        enrichments = []
        for group, df in res.items():
            df.insert(loc = 0, column = 'cell_classification', value = group)
            enrichments.append(df)
        enrichments = pd.concat(enrichments)
        enrichments.to_csv(out_enrichment, index = True, sep = '\t')
    
    return enrichments

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser('This script runs cell-type enrichments using scDRS')
    parser.add_argument('-i','--in', dest = '_in', help = 'Input gene list and weights', required = True)
    parser.add_argument('--h5ad', help = 'Input h5ad single-cell multiomics dataset', required = True)
    parser.add_argument('--label', nargs = '*', help = 'Columns containing cell classifications/types in the h5ad dataset',
        default = ['ROIGroup', 'ROIGroupCoarse', 'ROIGroupFine', 'roi', 'supercluster_term', 'cluster_id', 'subcluster_id', 'development_stage'])
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