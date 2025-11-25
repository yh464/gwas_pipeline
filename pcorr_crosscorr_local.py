#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-01-09

Creates a plot of phenotypic correlations at local level

Requires following inputs: 
    Phenotype file ready for GWAS analysis (with FID and IID)
'''

def corresp_global(loc, glob):
    from fnmatch import fnmatch
    if fnmatch(glob, '*_l_*') or fnmatch(glob, '*_r_*') or \
        fnmatch(glob, '*_l') or fnmatch(glob, '*_r'): return False
    glob = glob.replace('global','').replace('meta','')
    loc = loc.replace('local','').replace('regional','')
    return glob.find(loc) > -1
    
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
    from ._utils.path import normaliser
    norm = normaliser()
    
    for g2 in args.p2:
        print(f'Processing phenotype: {g2}')
        prefix = f'{args.out}/pcorr_local_{g2}.'+'_'.join(args.p1)
        
        summary = []
        df2 = pd.read_table(f'{args._in}/{g2}.txt', index_col = ['FID','IID'])
        for g1 in args.p1:
            df1 = pd.read_table(f'{args._in}/{g1}.txt', sep = '\\s+', index_col = ['FID','IID'])
            df_m = pd.concat([df1, df2], axis = 1)
            for p1 in df1.columns:
                for p2 in df2.columns:
                    if not compatible_roi(p1, p2) and not corresp_global(g1, p2): continue
                    tmp = df_m[[p1,p2]].dropna()
                    tmp.columns = [g1, g2]
                    r = tmp[g1].corr(tmp[g2])
                    n = tmp.shape[0]
                    se = ((1-r**2)/(n-2)) **0.5
                    z = r/se
                    p = sts.norm.cdf(-abs(z)) * 2
                    summary.append(
                        pd.DataFrame(dict(
                            group1 = g1, pheno1 = p1, group2 = g2, pheno2 = p2,
                            r = r, p = [p]))
                        )
        summary = pd.concat(summary)
        norm.normalise(summary).to_csv(f'{prefix}.txt', sep = '\t', index = False)
        wide = summary.pivot_table(index = 'pheno1', columns = 'group1', values = 'r')
        wide.columns.name = None
        wide.index.name = 'label'
        norm.normalise(wide).to_csv(f'{prefix}.wide.txt', sep = '\t', header = True, index = True)
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'This programme parses genetic correlation outputs and generates heatmaps')
    parser.add_argument('-p1', help = 'Group of local phenotypes', nargs = '*',
        default = ['deg_local','degi_local','degc_local', 'eff_local', 'clu_local','mpl_local'])
    parser.add_argument('-p2', help = 'Correlate phenotype, local or global', nargs = '*')
    parser.add_argument('-i','--in', dest = '_in', help = 'Phenotype file directory',
        default = '../pheno/ukb/')
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
        default = '../local_corr/')
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
        default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['_in','out']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    
    from ._utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/*',__file__)
    proj.add_output(args.out+'/pcorr_local_%pheno..*', __file__) # .* is a wildcard
    try: main(args)
    except: cmdhistory.errlog()