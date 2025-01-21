#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-01-08

Creates a plot of phenotypic correlations

Requires following inputs: 
    Phenotype file ready for GWAS analysis (with FID and IID)
'''

def main(args):
    import pandas as pd
    from _plots import corr_heatmap
    import scipy.stats as sts
    
    prefix = f'{args.out}/pcorr_'+'_'.join(args.p1)+'.'+ '_'.join(args.p2) 
    
    summary = []
    for g1 in args.p1:
        df1 = pd.read_table(f'{args._in}/{g1}.txt', sep = '\s+')
        for g2 in args.p2:
            df2 = pd.read_table(f'{args._in}/{g2}.txt', sep = '\s+')
            merge = pd.merge(df1, df2, on = ['FID','IID'])
            for t1 in df1.columns[2:]: # excludes FID and IID
                for t2 in df2.columns[2:]:
                    tmp = merge[[t1, t2]].dropna()
                    r = tmp[t1].corr(tmp[t2])
                    n = tmp.shape[0]
                    se = ((1-r**2)/(n-2)) **0.5
                    z = r/se
                    p = sts.norm.cdf(-abs(z)) * 2
                    summary.append(
                        pd.DataFrame(dict(
                            group1 = g1, pheno1 = t1,
                            group2 = g2, pheno2 = t2,
                            r = r, p = [p]))
                        )
    
    summary = pd.concat(summary)
    summary.to_csv(f'{prefix}.txt', sep = '\t', index = False)
    summary.pivot(index = ['group1','pheno1'], columns = ['group2','pheno2'], values = 'r'
        ).to_csv(f'{prefix}.wide.txt', index_label = False, sep = '\t', header = True, index = True)
    fig = corr_heatmap(summary)
    fig.savefig(f'{prefix}.pdf',
                bbox_inches = 'tight')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'This programme parses genetic correlation outputs and generates heatmaps')
    parser.add_argument('-p1', help = 'First group of phenotypes (x axis)', nargs = '*')
    parser.add_argument('-p2', help = 'Second group of phenotypes (y axis)', nargs = '*')
    parser.add_argument('-i','--in', dest = '_in', help = 'Phenotype file directory',
      default = '../pheno/ukb/')
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
      default = '../gene_corr/')
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
      default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['_in','out']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/*',__file__)
    proj.add_output(args.out+'/pcorr_%pheno..*', __file__) # .* is a wildcard
    try: main(args)
    except: cmdhistory.errlog()