#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-01-09

Creates a plot of phenotypic correlations at local level

Requires following inputs: 
    Phenotype file ready for GWAS analysis (with FID and IID)
'''

def compatible_roi(a, b):
    a = a.lower(); b = b.lower()
    if a[0] == 'x': a = a[1:] # unilateral regions starting with a number may be masked with 'x'
    if b[0] == 'x': b = b[1:]
    a = a.replace('.','_').replace('-','_')
    if a in ['fid', 'iid', 'eid'] or b in ['fid', 'iid', 'eid']:
        return False
    a = a.replace('_roi','').replace('_0.01',''); b = b.replace('_roi','').replace('_0.01','')
    if a[:2] in ['lh','l_','rh','r_'] and b[:2] in ['lh','l_','rh','r_']: # unilateral regions
        return a == b
    
    # if not both 
    def unil(a):
        if a[:3] in ['lh_','rh_']: a = a[3:]
        if a[:2] in ['l_','r_']: a = a[2:]
        return a
    
    return unil(a) == unil(b)

def main(args):
    import pandas as pd
    import scipy.stats as sts
    
    for g2 in args.p2:
        print(f'Processing phenotype: {g2}')
        prefix = f'{args.out}/pcorr_local_{g2}.'+'_'.join(args.p1)
        
        summary = []
        df2 = pd.read_table(f'{args._in}/{g2}.txt', sep = '\s+')
        for g1 in args.p1:
            df1 = pd.read_table(f'{args._in}/{g1}.txt', sep = '\s+')
            eid = pd.merge(df1[['FID','IID']], df2[['FID','IID']])
            df1m = pd.merge(df1, eid); df2m = pd.merge(df2, eid) # this design is to prevent duplicate column names
            for t1 in df1.columns[2:]: # excludes FID and IID
                for t2 in df2.columns[2:]:
                    if not compatible_roi(t1, t2): continue
                    tmp = pd.concat([df1m[t1], df2m[t2]],axis = 1).dropna()
                    tmp.columns = [g1, g2]
                    r = tmp[g1].corr(tmp[g2])
                    n = tmp.shape[0]
                    se = ((1-r**2)/(n-2)) **0.5
                    z = r/se
                    p = sts.norm.cdf(-abs(z)) * 2
                    summary.append(
                        pd.DataFrame(dict(
                            pheno = g1, roi = t1,
                            r = r, p = [p]))
                        )
        summary = pd.concat(summary)
        summary.to_csv(f'{prefix}.txt', sep = '\t', index = False)
        wide = summary.pivot(index = 'roi', columns = 'pheno', values = 'r')
        wide.columns.name = None
        wide.index.name = 'label'
        wide.to_csv(f'{prefix}.wide.txt', sep = '\t', header = True, index = True)
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'This programme parses genetic correlation outputs and generates heatmaps')
    parser.add_argument('-p1', help = 'Group of local phenotypes', nargs = '*',
        default = ['deg_local','degi_local','degc_local', 'eff_local', 'clu_local','mpl_local'])
    parser.add_argument('-p2', help = 'Correlate local phenotype', nargs = '*')
    parser.add_argument('-i','--in', dest = '_in', help = 'Phenotype file directory',
        default = '../pheno/ukb/')
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
        default = '../local_corr/')
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
    proj.add_output(args.out+'/pcorr_local_%pheno..*', __file__) # .* is a wildcard
    try: main(args)
    except: cmdhistory.errlog()