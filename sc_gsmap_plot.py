#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2025-09-22

Plots spatial enrichment results from gsMap in a scatterplot

Requires following inputs: 
    gsMap spatial LDSC output files (.csv.gz)
    original spatial transcriptomics dataset for gsMap (.h5ad)
'''

def main(args):
    import os
    from _utils.path import find_gwas
    from _plots import colourcode_scatterplot
    from _plots.aes import redgrey
    import matplotlib.pyplot as plt
    import scanpy as sc
    import pandas as pd
    import warnings
    import numpy as np
    from tqdm import tqdm

    # find ST datasets and phenotype files
    st_datasets = os.listdir(args.st)
    pheno = find_gwas(args.pheno, dirname = args.gwa, ext = 'sumstats', long = True)

    for g, p in pheno:
      for s in tqdm(st_datasets, desc = f'{g}/{p}'):
        h5ad = f'{args.st}/{s}/find_latent_representations/{s}_add_latent.h5ad'
        gsmap_output = f'{args._in}/{g}/{p}/{s}_spatial_ldsc.csv.gz'
        out_fig = f'{args._in}/{g}/{p}/{s}_gsmap_spatial_ldsc.png'

        if not os.path.isfile(h5ad): warnings.warn(f'Missing h5ad file for {s}, skipping'); continue
        if not os.path.isfile(gsmap_output):
          warnings.warn(f'Missing spatial LDSC output file, please run:\n    python sc_gsmap_batch.py {g}/{p}')
          continue
        if os.path.isfile(out_fig) and not args.force: continue

        adata = sc.read_h5ad(h5ad, 'r')
        coords = pd.DataFrame(index = adata.obs_names, columns = ['x','y'], data = adata.obsm['spatial'])
        sldsc = pd.read_table(gsmap_output, sep = ',', index_col = 'spot', dtype = {'spot':str})
        df = coords.join(sldsc, how = 'inner')
        if df.shape[0] == 0: warnings.warn(f'No overlapping spots between {h5ad} and {gsmap_output}, skipping'); continue
        df['logp'] = -np.log10(df['p'])
        colourcode_scatterplot.scatterplot_noaxis(df['x'], df['y'], df['logp'], redgrey, 0.1, rep = False, vname = r"$-log_{10}{(P)}$")
        plt.savefig(out_fig, dpi = 400)
        plt.close()

if __name__ == '__main__':  
    from argparse import ArgumentParser
    parser = ArgumentParser(description = 'This script runs cell-type enrichments using scDRS')
    parser.add_argument('pheno', nargs = '*', help = 'Phenotypes')
    parser.add_argument('-i', '--in', dest = '_in', help = 'input gsMap result directory', default = '../sc/gsmap')
    parser.add_argument('-g','--gwa', help = 'Directory containing LDSC summary statistics',
        default = '../gcorr/ldsc_sumstats')
    parser.add_argument('-s','--st', help = 'Directory containing gsMap processed spatial transcriptomics data',
        default = '/rds/project/rds-Nl99R8pHODQ/multiomics/gsmap') # intentionally absolute
    parser.add_argument('-f','--force', help = 'force overwrite', default = False, action = 'store_true')
    args = parser.parse_args()
    
    # path normalisation
    import os
    for arg in ['gwa','st','_in']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))

    from _utils import logger, cmdhistory
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()