#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-03-24

Batch submits jobs for MVMR for all GWAS files in a directory. 
Scans the entire directory for GWAS summary stats of the same data extension.
'''

def main(args):
    import pandas as pd
    import scipy.stats as sts
    import os
    from logparser import parse_rg_log
    from _utils.path import find_gwas
    
    force = '-f' if args.force else ''
    # array submitter
    from _utils import array_submitter
    submitter = array_submitter.array_submitter(
        name = 'mvmr_'+'_'.join(args.p2), env = 'gentoolsr',
        n_cpu = 4, 
        timeout = 60,
        debug = True
        )
    
    # screen for GWA summary stats
    exposures = find_gwas(args.p1, args.gwa, args.e1)
    outcomes = find_gwas(args.p2, args.gwa, args.e2)
    
    # for each outcome
    for p2,x2 in outcomes:
        # create output directory
        outdir = f'{args.out}/{p2}/mvmr'
        if not os.path.isdir(outdir):os.system(f'mkdir -p {outdir}')
        out_prefix = f'{args.out}/{p2}/mvmr/{x2}_mvmr_' + ('_'.join(args.p1)).replace('_local','')
        if os.path.isfile(f'{out_prefix}.txt') and not args.force: continue
        
        # screen for only correlated exposures
        exposures_corr = []
        for p1, x1 in exposures:
            if p1 < p2: fname = f'{args.rg}/{p1}.{p2}/{p1}_{x1}.{p2}_{x2}.rg.log'
            else: fname = f'{args.rg}/{p2}.{p1}/{p2}_{x2}.{p1}_{x1}.rg.log'
            
            if not os.path.isfile(fname): continue
            rg, se = parse_rg_log(fname)
            p = 1-sts.chi2.cdf((rg/se)**2, df = 1) # p value
            exposures_corr.append(pd.DataFrame(group = p1, pheno = x1, pval= [p]))
        exposures_corr = pd.concat(exposures_corr)
        exposures_corr['q'] = 1
        for p1,_ in exposures:
            exposures_corr.loc[(exposures_corr.group==p1) & (~exposures.corr.pval.isna()),'q'] = \
                sts.false_discovery_control(exposures_corr.loc[
                (exposures_corr.group==p1) & (~exposures.corr.pval.isna()),'pval'])
        exposures_corr = exposures_corr.loc[exposures_corr.q < 0.05,['group','pheno']]
        exposures_corr['exp'] = exposures_corr['group'] + '/' + exposures_corr['pheno']
        exposures_filtered = exposures_corr.exp.to_list()
        
        # find instruments
        instruments = []
        exposures_corr = exposures_corr['group'].unique()
        for i in range(len(exposures_corr)):
            p1i = exposures_corr[i]
            for j in range(i, len(exposures_corr)):
                p1j = exposures_corr[j]
                instruments.append(f'{args._in}/{p1j}_clumped_for_{p1i}_{args.pval:.0e}.txt')
            instruments.append(f'{args._in}_{p2}_clumped_for_{p1i}_{args.pval:.0e}.txt')
        
        cmd = ['Rscript mr_mvmr.r -i'] + instruments + ['--p1'] + exposures_filtered + \
            [f'--p2 {p2}/{x2} --pval {args.pval} -o {out_prefix}',force]
        submitter.add(' '.join(cmd))
    submitter.submit()
    return submitter

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'This script batch runs MR for groups of phenotypes')
    path_spec = parser.add_argument_group('Path specifications')
    path_spec.add_argument('-g','--gwa', dest = 'gwa', 
                        help = 'GWA directory, assumes both groups of pheno to be in the same dir',
                        default = '../gwa')
    path_spec.add_argument('-i','--in', dest = '_in', 
                        help = 'Directory of instruments, output of mr_extract_snp_batch.py',
                        default = '../mr/instruments')
    path_spec.add_argument('-c','--clump', dest = 'clump', help = 'Directory of clumping files',
                        default = '../clump')
    path_spec.add_argument('-rg', dest = 'rg', help = 'Directory to rg log files',
                        default = '../gcorr/rglog')
    path_spec.add_argument('-o','--out', dest = 'out', help = 'Output directory',
                        default = '../mr')
    pheno_spec = parser.add_argument_group('Phenotype specifications')
    pheno_spec.add_argument('-p1','--pheno1', dest = 'p1', 
                        help = 'Exposure', nargs = '*', default = [])
    pheno_spec.add_argument('-e1','--ext1', dest = 'ext1', help = 'Extension for phenotype group 1',
                        default = 'fastGWA')
    pheno_spec.add_argument('-p2','--pheno2', dest = 'p2', 
                        help = 'Outcome', nargs = '*', 
                        default = ['disorders', 'disorders_subtypes'])
    pheno_spec.add_argument('-e2','--ext2', dest = 'ext2', help = 'Extension for phenotype group 2',
                        default = 'fastGWA')
    parser.add_argument('--pval', help = 'Clumping p-value threshold', default = 5e-8, type = float)
    parser.add_argument('-f','--force', dest = 'force', action = 'store_true',
                        default = False, help = 'Force overwrite')
    args = parser.parse_args()
    args.p1.sort(); args.p2.sort()
    
    import os
    for arg in ['_in','out','clump', 'h2','rg']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    try: args.n1 = int(args.n1)
    except: args.n1 = os.path.realpath(args.n1)
    try: args.n2 = int(args.n2)
    except: args.n2 = os.path.realpath(args.n2)
    
    from _utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_input(f'{args._in}/{args.p1}/*.{args.ext1}', __file__)
    proj.add_input(f'{args._in}/{args.p2}/*.{args.ext2}', __file__)
    proj.add_input(f'{args.clump}/{args.p1}/*.clumped',__file__)
    proj.add_input(f'{args.clump}/{args.p2}/*.clumped',__file__)
    proj.add_output(f'{args.out}/{args.p2}/*',__file__)
    
    try: main(args)
    except: cmdhistory.errlog()