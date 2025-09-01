#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2025-09-01

A python implementation of cepo to save computing time

Requires following inputs: 
    h5ad file and columns for cell classifications
'''

from ast import expr
from html import parser
from math import sin


def segidx(mat):
    import pandas as pd
    import numpy as np
    nz = pd.Series((mat > 0).mean(axis = 1))
    ms = mat.mean(axis = 1)
    sds = np.std(mat, axis = 1)
    cvs = pd.Series(sds/ms)
    
    x1 = nz.rank()/(len(nz)+1)
    x2 = 1 - cvs.rank()/(len(cvs)+1)

    return (x1 + x2)/2

def singlebatchcepo(expr, genes, celltypes, mincells = 20, exprspct:float = 0.05):
    # columns = cells, rows = genes
    import pandas as pd
    celltypes = pd.Series(celltypes).astype('category')
    celltypes_df = pd.get_dummies(celltypes); del celltypes
    for ct in celltypes_df.columns:
        if celltypes_df[ct].sum() < mincells:
            celltypes_df = celltypes_df.drop(columns = ct)
    if celltypes_df.shape[1] == 0: return []

    # filter genes such that each gene is expressed in at least exprspct of cells in at least one cell type
    keep = pd.DataFrame(index = genes, columns = celltypes_df.columns, data = False)
    for ct in celltypes_df.columns:
        keep[ct] = (expr[:, celltypes_df[ct]] > 0).sum(axis = 1) >= celltypes_df[ct].sum() * exprspct
    keep = keep.any(axis = 1)

    segmat = pd.DataFrame(index = genes[keep], columns = celltypes_df.columns, data = float(0))
    for ct in celltypes_df.columns:
        segmat.loc[keep, ct] = segidx(expr[keep, celltypes_df[ct]])
    
    out = pd.DataFrame(index = genes, columns = celltypes_df.columns, data = float(0))
    for ct in celltypes_df.columns:
        out.loc[keep,ct] = (segmat[ct] - segmat.drop(ct,axis=1)).mean(axis = 1)
    return out


def main(args):
    import scanpy as sc
    import pandas as pd
    import scipy.sparse as sp
    import os

    if os.path.isfile(f'{args.out}.cepo.txt') and not args.force: return

    adata = sc.read_h5ad(args._in)
    adata.X = sp.csc_array(adata.X)
    genes = adata.var.index
    cells = adata.obs.index

    out_dfs = []
    for label in args.label:
        if label not in adata.obs.columns: Warning(f'Label {label} not found in adata.obs, skipping'); continue

        df = singlebatchcepo(adata.X, genes, adata.obs[label], mincells = 20, exprspct = 0.05)
        df.columns = [f'cepo.{label}.{x}' for x in df.columns]
        out_dfs.append(df)
        del df
    out_dfs = pd.concat(out_dfs, axis = 1)
    out_dfs.index.name = 'gene'
    out_dfs.to_csv(f'{args.out}.cepo.txt', sep = '\t', index = True, header = True)
    return out_dfs

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description = 'This script generates gene scores for MAGMA gene covariate analysis')
    parser.add_argument('-i','--in', dest = '_in', help = 'Input h5ad single-cell multiomics dataset', required = True)
    parser.add_argument('-o','--out', help = 'Output prefix', required = True)
    parser.add_argument('--label', nargs = '*', help = 'Columns containing cell classifications/annotations in the h5ad dataset',
        default = ['ROIGroup','ROIGroupCoarse', 'ROIGroupFine', 'roi', 'supercluster_term', 'cluster_id', 'subcluster_id', 'development_stage'])
    parser.add_argument('-f','--force', help = 'Force overwrite', action = 'store_true', default = False)
    args = parser.parse_args()

    import os
    args._in = os.path.realpath(args._in); args.out = os.path.realpath(args.out)

    from _utils import logger, cmdhistory
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()