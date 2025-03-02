#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2023-07-13

Creates a plot of genetic, phenotypic correlations and heritability

Requires following inputs: 
    Phenotype files
    LDSC rg logs
    LDSC h2 logs
'''

def main(args):
    import pandas as pd
    import numpy as np
    import scipy.stats as sts
    from _utils.path import normaliser
    from _plots import corr_heatmap
    
    os.chdir(args._in)
    
    # Phenotypic matrix
    pheno_df = []
    pheno_list = []
    for x in args.pheno:
        df = pd.read_table(f'{x}.txt', sep = '\\s+', index_col = ['FID','IID'])
        pheno_list += [x] * df.shape[1]
        pheno_df.append(df)
    pheno_df = pd.concat(pheno_df, axis = 1)
    prefix_list = pheno_df.columns.tolist()
    del df
    n = pheno_df.shape[0]
    
    # Correlation matrices
    pcorr_summary = []
    gcorr_summary = []
    for x in range(len(prefix_list)):
        p1 = pheno_list[x]; x1 = prefix_list[x]
        if any([x1.find(ex) > -1 for ex in args.exclude]): continue
    
        fname = f'{args.h2}/{p1}/{x1}.h2.log'
        if os.path.isfile(fname):
            tmp = open(fname).read().splitlines()[-7].replace('(','').replace(')','').split()
            while tmp.count('') > 0:
                tmp.remove('')
            try: h2 = float(tmp[-2])
            except: h2 = np.nan
            try: se = max([float(tmp[-1]),1e-20])
            except: se = np.nan
            
            p = 1-sts.chi2.cdf((h2/se)**2, df = 1)
            gcorr_summary.append(pd.DataFrame(dict(group1 = p1, pheno1 = x1,
                group2 = p1, pheno2 = x1, r = h2, se = se, p = [p])))
            
        for y in range(x, len(prefix_list)):
            p2 = pheno_list[y]; x2 = prefix_list[y]
            if any([x2.find(ex) > -1 for ex in args.exclude]): continue
            
            if p1 < p2 or x1 < x2: tmp1, tmx1, tmp2, tmx2 = p1, x1, p2, x2
            else: tmp1, tmx1, tmp2, tmx2 = p2, x2, p1, x1
            
            # phenotypic correlation
            rp = pheno_df[x1].corr(pheno_df[x2]) # automatically drops NaN
            se = max([((1-rp**2)/(n-2))**0.5, 1e-20])
            pcorr_summary.append(pd.DataFrame(dict(
                group1 = [tmp1], pheno1 = tmx1, 
                group2 = tmp2, pheno2 = tmx2,
                r = rp, se = se, p = 1-sts.chi2.cdf((rp/se)**2, df = 1)
                )))
            
            # genetic correlation
            fname = f'{args.corr}/{tmp1}.{tmp2}/{tmp1}_{tmx1}.{tmp2}_{tmx2}.rg.log'
            if os.path.isfile(fname):
                tmp = open(fname)
                tmp_stats = tmp.read().splitlines()
                tmp_stats = tmp_stats[-4].split()
                while tmp_stats.count('') > 0:
                  tmp_stats.remove('')
                try: 
                    rg = float(tmp_stats[2])
                    if rg > 1: rg = 1
                    if rg < -1: rg = -1
                except: 
                  rg = np.nan
                  print(f'{fname} shows NA correlation!')
                try: se = max((float(tmp_stats[3]),1e-20))
                except: se = np.nan
            
                p = 1-sts.chi2.cdf((rg/se)**2, df = 1) # p value
                
                gcorr_summary.append(pd.DataFrame(dict(group1 = tmp2, pheno1 = tmx2,
                  group2 = tmp1, pheno2 = tmx1, r = rg, se = se, p = [p])))
            
    pcorr_summary = pd.concat(pcorr_summary)
    pcorr_summary['q'] = np.nan
    pcorr_summary.loc[~pcorr_summary.p.isna(),'q'] = sts.false_discovery_control(
        pcorr_summary.loc[~pcorr_summary.p.isna(),'p'])
    
    gcorr_summary = pd.concat(gcorr_summary)
    gcorr_summary['q'] = np.nan
    gcorr_summary.loc[~gcorr_summary.p.isna(),'q'] = sts.false_discovery_control(
        gcorr_summary.loc[~gcorr_summary.p.isna(),'p'])
    
    # plot heatmap
    prefix = '_'.join(args.pheno)
    norm = normaliser()
    summary = norm.normalise(pd.concat([pcorr_summary, gcorr_summary]))
    summary.to_csv(f'{args.out}/correlation_{prefix}.txt', sep = '\t', index = False)
    
    fig = corr_heatmap(summary, absmax = 1, autocor = True, annot = 'Heritability')
    # fig.text(0.5, 0, 'Phenotypic Correlation', size = 'x-large')
    # fig.text(0, 0.5, 'Genetic Correlation', size = 'x-large', rotation = 'vertical')
    fig.savefig(f'{args.out}/correlation_{prefix}.pdf', bbox_inches = 'tight')
    fig.savefig(f'{args.out}/correlation_{prefix}.png', bbox_inches = 'tight')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'This programme parses genetic correlation outputs and generates heatmaps')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
    parser.add_argument('-i','--in', dest = '_in', help = 'Phenotype file directory',
      default = '../pheno/ukb/')
    parser.add_argument('--corr', dest = 'corr', help = 'correlational logfile directory',
      default = '../gcorr/rglog/')
    parser.add_argument('--h2', dest = 'h2', help = 'Heritability logfile directory',
      default = '../gcorr/ldsc_sumstats/')
    parser.add_argument('--exclude', help = 'phenotypes to exclude', nargs = '*', default = [])
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
      default = '../gcorr/')
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
      default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','corr','h2']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    args.pheno.sort()
    
    from _utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args.corr+'/gcorr/%pheno_%maf.%pheno_%maf.rg.log', __file__)
    proj.add_input(args.h2+'/%pheng/%pheno_%maf.h2.log',__file__)
    proj.add_output(args.out+'/correlation_%pheno..*', __file__) # .* is a wildcard
    proj.add_output(args.out+'/h2_%pheno..*', __file__) # .* is a wildcard)
    try: main(args)
    except: cmdhistory.errlog()