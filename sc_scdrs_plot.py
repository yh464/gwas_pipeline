#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2025-08-27

Plots scDRS results in a heatmap

Requires following inputs: 
    scDRS output folder
    columns containing cell classifications
'''
def main(args):
    # find phenotypes
    from _utils.path import find_gwas
    pheno = find_gwas(args.pheno, long = True)
    pheno_short = find_gwas(args.pheno); 

    # find h5ad annotations
    import os
    for sc in args.sc:
        h5ad_prefix = [x[:-5] for x in os.listdir(f'{args.h5ad}/{sc}') if x[-5:] =='.h5ad']
        out_prefix = f'{args._in}/'+'_'.join([x[0] for x in pheno_short])+f'.{sc}'
        # concatenate scDRS output
        import pandas as pd
        summary = []
        for g, p in pheno:
            pheno_summary = []
            for h5prefix in h5ad_prefix:
                try: df = pd.read_table(f'{args._in}/{g}/{p}/{sc}/{p}.{h5prefix}.scdrs.enrichment.txt').assign(dataset = h5prefix)
                except: Warning(f'Missing scDRS enrichment for {p}.{sc}.{h5prefix}'); continue
                pheno_summary.append(df)
            if len(pheno_summary) == 0: Warning(f'Missing scDRS enrichment for {g}/{p}'); continue
            pheno_summary = pd.concat(pheno_summary, axis = 0)
            pheno_summary.to_csv(f'{args._in}/{g}/{p}.{sc}.scdrs.enrichment.txt', index = False, sep = '\t')
            pheno_summary.columns = ['cell_type', 'annotation', 'n_cell','n_ctrl','p','beta','hetero_p','hetero_z','fdr.05','fdr.1','fdr.2','dataset']
            pheno_summary = pheno_summary.assign(group = g, pheno = p)
            pheno_summary = pheno_summary[['group','pheno','dataset','cell_type','beta','p', # necessary columns for _plots.corr_heatmap
                'fdr.05','fdr.1','fdr.2','annotation','n_cell','n_ctrl','hetero_p','hetero_z']]
            summary.append(pheno_summary)
        
        # plot heatmap
        summary = pd.concat(summary)
        from _plots.corr_heatmap import corr_heatmap
        for lab in args.label:
            summary.loc[summary.annotation == lab,:].to_csv(f'{out_prefix}_{lab}_enrichment.txt', index = False, sep = '\t')
            fig = corr_heatmap(summary.loc[summary.annotation == lab,:], p_threshold = [0.05, 0.001])
            fig.savefig(f'{out_prefix}_{lab}_enrichment.pdf', bbox_inches = 'tight')


if __name__ == '__main__':  
    from argparse import ArgumentParser
    parser = ArgumentParser(description = 'This script runs cell-type enrichments using scDRS')
    parser.add_argument('pheno', nargs = '*', help = 'Phenotypes')
    parser.add_argument('-i','--in', dest = '_in', help = 'Directory containing scDRS output', default = '../sc/scdrs')
    parser.add_argument('-s','--sc', nargs = '*', help = 'single-cell dataset', default = ['siletti_2023','wang_2025'])
    parser.add_argument('--h5ad', help = 'Input directory containing h5ad single-cell multiomics dataset',
        default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/multiomics/scdrs') # intentionally absolute
    parser.add_argument('--label', nargs = '*', help = 'Columns containing cell classifications/annotations in the h5ad dataset',
        # default = ['ROIGroup', 'ROIGroupCoarse', 'ROIGroupFine', 'roi', 'supercluster_term', 'cluster_id', 'subcluster_id', 'development_stage'])
        default = ['supercluster_term', 'cluster_id', 'subcluster_id'])
    parser.add_argument('-f','--force',dest = 'force', help = 'force overwrite', default = False, action = 'store_true')
    args = parser.parse_args()
    
    # path normalisation
    import os
    for arg in ['_in','h5ad']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))

    from _utils import logger, cmdhistory
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()