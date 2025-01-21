#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-01-10

Creates a plot of PRS correlations based on existing summary files

Required workflow:
    prs_corr_batch.py
'''

def main(args):
    import os
    import pandas as pd
    os.chdir(args._in)
    
    out_fname = f'{args.out}/prscorr_'+'_'.join(args.pheno)
    
    # read all summary files
    summary = []
    for p in args.pheno:
        f = f'{args._in}/{p}_summary.txt'
        print(f'Reading {f}')
        summary.append(pd.read_table(f))
    summary = pd.concat(summary)
    
    from _plots import corr_heatmap
    fig = corr_heatmap(summary)
    fig.savefig(f'{out_fname}.pdf', bbox_inches = 'tight')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Computes correlational plots between IDPs and PGS (batch)')
    parser.add_argument('pheno', nargs = '*', help = 'Phenotype files to process')
    parser.add_argument('-i','--in', dest = '_in', help = 'Input directory',
      default = '../prs/prs_corr/')
    parser.add_argument('-o','--out', dest = 'out', help = 'Output directory')
    # always overwrites
    args=parser.parse_args()
    import os
    args._in = os.path.realpath(args._in)
    if args.out == None: args.out = args._in
    else: args.out = os.path.realpath(args.out)
        
    from _utils import cmdhistory, logger
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()