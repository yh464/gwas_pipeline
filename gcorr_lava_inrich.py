#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-11-02

Conducts down-stream interval-based enrichment for LAVA results

Requires following inputs: 
    LAVA output
    inrich default resources
Outputs:
    Main analysis (pathway enrichments)
    Interval-gene-target analysis
'''

import os

def parse_lava_output(file):
    import pandas as pd
    df = pd.read_table(file).rename(columns = {'CHR':'chr', 'START': 'start_bp', 'STOP':'stop_bp'})
    df = df.loc[df['type'] == 'rg', :].reset_index(drop = True)
    out = []
    for phen, phen_df in df.groupby(['group1','pheno1','group2','pheno2']):
        phen_dict = {'+': dict(), '-': dict()}
        phen_df = phen_df.loc[phen_df['p'] < 0.05, :].reset_index(drop = True)
        if any(phen_df.stat > 0): phen_dict['+']['p0.05'] = phen_df.loc[phen_df['stat'] > 0, ['chr','start_bp','stop_bp']]
        if any(phen_df.stat < 0): phen_dict['-']['p0.05'] = phen_df.loc[phen_df['stat'] < 0, ['chr','start_bp','stop_bp']]
        phen_df = phen_df.loc[phen_df['q'] < 0.1, :].reset_index(drop = True)
        if any(phen_df.stat > 0): phen_dict['+']['fdr0.1'] = phen_df.loc[phen_df['stat'] > 0, ['chr','start_bp','stop_bp']]
        if any(phen_df.stat < 0): phen_dict['-']['fdr0.1'] = phen_df.loc[phen_df['stat'] < 0, ['chr','start_bp','stop_bp']]
        phen_df = phen_df.loc[phen_df['q'] < 0.05, :].reset_index(drop = True)
        if any(phen_df.stat > 0): phen_dict['+']['fdr0.05'] = phen_df.loc[phen_df['stat'] > 0, ['chr','start_bp','stop_bp']]
        if any(phen_df.stat < 0): phen_dict['-']['fdr0.05'] = phen_df.loc[phen_df['stat'] < 0, ['chr','start_bp','stop_bp']]
        out.append((phen[0], phen[1], phen[2], phen[3], phen_dict))
    return out

def main(args):
    import pandas as pd
    from _utils.plugins.inrich import inrich
    lava_results = parse_lava_output(args._in)
    all_main = []; all_igt = []

    for g1, p1, g2, p2, phen_dict in lava_results:
        for sign in ['+','-']:
            for threshold, df in phen_dict[sign].items():
                if df.shape[0] == 0: continue
                main_analysis, igt_analysis = inrich(df, chrom_col='chr', start_col='start_bp', stop_col='stop_bp',
                    inrich_dir = args.inrich)
                if main_analysis.shape[0] > 0:
                    main_analysis = main_analysis.assign(
                        group1 = g1, pheno1 = p1, group2 = g2, pheno2 = p2,
                        sign = sign, threshold = threshold)
                    all_main.append(main_analysis)
                if igt_analysis.shape[0] > 0:
                    igt_analysis = igt_analysis.assign(
                        group1 = g1, pheno1 = p1, group2 = g2, pheno2 = p2,
                        sign = sign, threshold = threshold)
                    all_igt.append(igt_analysis)
    if len(all_main) > 0:
        all_main = pd.concat(all_main, axis = 0)
        all_main.to_csv(f'{args.out}.txt', sep = '\t', index = False)
    if len(all_igt) > 0:
        all_igt = pd.concat(all_igt, axis = 0)
        all_igt.to_csv(f'{args.out}.igt.txt', sep = '\t', index = False)
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'This scripts investigates interval-based enrichments for LAVA-identified genomic regions')
    parser.add_argument('-i', '--in', dest = '_in', type = str, required = True, help = 'LAVA output file')
    parser.add_argument('--inrich', 
        help = 'folder of the inrich binary and resources, should contain resources/genes.txt and resources/snps.txt',
        default = '/rds/project/rds-Nl99R8pHODQ/toolbox/inrich')
    parser.add_argument('-o', '--out', type = str, required = True, help = 'output prefix')
    args = parser.parse_args()
    for key, value in vars(args).items(): setattr(args, key, os.path.realpath(value))
    main(args)