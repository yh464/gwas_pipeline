#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2025-10-14

Conducts scDRS scoring for a single phenotype

Requires following inputs: 
    Gene-level summary statistics
        Column specifications: gene ID (ENSGxxxxxxxx), p-value and effect size (optional)
    P-value threshold for sig genes (default fdr < 0.05)
    List of number of genes to test as gene sets
    Pre-processed scDRS single cell data
        Column specifications: Columns of adata.obs that contain cell group labels
            (for Siletti et al 2023: ROIGroup, ROIGroupCoarse, ROIGroupFine, roi, supercluster_term, cluster_id, subcluster_id, development_stage)
'''

def main(args = None, **kwargs):
    from _utils.gadgets import namespace
    import os
    import scanpy as sc
    import scdrs
    import pandas as pd
    from tqdm import tqdm
    if args == None:
        from _utils.gadgets import namespace
        args = namespace(**kwargs)

    out_score = f'{args.out}.sensitivity.txt'
    if os.path.isfile(out_score) and not args.force: return
    weights = pd.read_table(args._in).sort_values('weight', ascending = False).reset_index(drop = True)
    gene_list = weights['gene'].iloc[:max(args.nsig)].tolist()
    gene_weight = weights['weight'].iloc[:max(args.nsig)].values if 'weight' in weights.columns else None
    adata = sc.read_h5ad(args.h5ad, 'r')
    block_size = int(2e9/adata.shape[1])
    if adata.shape[0] >= block_size:
        score_df = []
        for i in range(0, adata.shape[0], block_size):
            # control genes are the same as the same random seed is used
            cmin = i; cmax = min(i+block_size, adata.shape[0])
            tempfile = f'{args.out}.sensitivity.{cmin}_{cmax}.tmp.gz'
            try: 
                if args.force: raise FileNotFoundError
                temp_score = pd.read_table(tempfile, index_col = 0)
            except:
                print(f'Scoring cells {cmin} - {cmax}')
                temp = adata[cmin:cmax,:]
                temp_score = pd.DataFrame(0., index = temp.obs_names, columns = [f'n{x}' for x in args.nsig])
                for nsig in tqdm(args.nsig, desc = f'Scoring with different number of genes'):
                    tmpdf = scdrs.score_cell(temp, gene_list[:nsig], 
                        gene_weight[:nsig] if gene_weight is not None else None, n_ctrl = 1,
                        return_ctrl_norm_score = False, verbose = False)
                    print(tmpdf)
                    temp_score[f'n{nsig}'] = tmpdf['raw_score']
                temp_score.to_csv(tempfile, index = True, sep = '\t')
            score_df.append(temp_score)
        score = pd.concat(score_df)
    else: score = scdrs.score_cell(adata, gene_list, gene_weight, return_ctrl_norm_score = True, verbose = True)
    score.to_csv(out_score, index = True, sep = '\t')
    score.corr().to_csv(f'{args.out}.sensitivity.corr.txt', index = True, sep = '\t')
    adata.file.close()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser('This script runs cell-type enrichments using scDRS')
    parser.add_argument('-i','--in', dest = '_in', help = 'Input gene list and weights', required = True)
    parser.add_argument('--h5ad', help = 'Input h5ad single-cell multiomics dataset', required = True)
    parser.add_argument('--label', nargs = '*', help = 'Columns containing cell classifications/annotations in the h5ad dataset',
        default = ['ROIGroup', 'ROIGroupCoarse', 'ROIGroupFine', 'roi', 'supercluster_term', 'cluster_id', 'subcluster_id', 'development_stage', # siletti
        'Class','Subclass','Type_updated', 'Cluster', 'Tissue', # wang
        'subcluster_identity_broad','subcluster_identity', # keefe
        ])
    parser.add_argument('-n','--nsig', help = 'Number of significant genes', nargs = '*', type = int,
        default = [100, 200, 500, 1000, 2000])
    parser.add_argument('-d', '--downstream', help = 'Conduct downstream analyses', default = False, action = 'store_true')
    parser.add_argument('-o', '--out', dest = 'out', help = 'output prefix, no .txt', required = True)
    parser.add_argument('-f','--force',dest = 'force', help = 'force overwrite', default = False, action = 'store_true')
    args = parser.parse_args()
    
    # path normalisation
    import os
    for arg in ['_in','h5ad','out']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    args.nsig = sorted(list(set(args.nsig))) # only unique elements

    from _utils import cmdhistory
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()