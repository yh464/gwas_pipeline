#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-03-24

Batch submits jobs for MVMR for all GWAS files in a directory. 
Scans the entire directory for GWAS summary stats of the same data extension.
'''

def main(args):
    from ._plugins.logparser import crosscorr_parse
    from ._utils.path import find_gwas
    import warnings
    
    # extract SNPs
    from mr_extract_snp_batch import api
    snp_submitter = api(p1 = args.p1, p2 = args.p2, _in = args.gwa,
                        out = args._in, clump = args.clump)
    
    force = '-f' if args.force else ''
    
    # screen for GWA summary stats
    exposures = find_gwas(args.p1, dirname=args.gwa, ext=args.ext1, se = True)
    outcomes = find_gwas(args.p2, dirname=args.gwa, ext=args.ext2, se = True)
    
    # array submitter
    timeout = min(sum([len(x[1]) for x in exposures])*5,720)
    from ._utils.slurm import array_submitter
    submitter = array_submitter(
        name = f'mvmr_{args.p1[0]}.'+'_'.join(args.p2), env = 'gentoolsr',
        n_cpu = 1, 
        timeout = timeout,
        dependency = snp_submitter,
        )
    
    # correlation matrix between exposures and outcomes
    exp_corr_out = crosscorr_parse(exposures, outcomes, logdir=args.rg)
    exp_corr_exp = crosscorr_parse(exposures, logdir = args.rg, h2dir = None)
    exp_corrmat = os.path.realpath(f'{args.rg}/../corr_'+'_'.join(args.p1) + '.txt')
    exp_corr_exp.to_csv(exp_corrmat, sep = '\t', index = False)
    
    # for each outcome
    for g2,p2 in [(x,y) for x,z in outcomes for y in z]:
        # create output directory
        outdir = f'{args.out}/{g2}/mvmr'
        if not os.path.isdir(outdir):os.system(f'mkdir -p {outdir}')
        out_prefix = f'{args.out}/{g2}/mvmr/{p2}_mvmr_' + ('_'.join(args.p1)).replace('_local','')
        all_output = ['_mvmrivw.txt', '_mvmrhorse.txt','_mvmrcmlsusie.txt']
        if all([os.path.isfile(f'{out_prefix}{z}') for z in all_output]) and \
            not args.force: continue
        
        # screen for only correlated exposures
        if args.rgp < 0:
            exposures_corr = exp_corr_out.loc[(exp_corr_out.group2==g2) & (exp_corr_out.pheno2==p2) &\
                (exp_corr_out.q < 0.05),['group1','pheno1']]
        else:
            exposures_corr = exp_corr_out.loc[(exp_corr_out.group2==g2) & (exp_corr_out.pheno2==p2) &\
                (exp_corr_out.p < args.rgp),['group1','pheno1']]
        exposures_corr['exp'] = exposures_corr['group1'] + '/' + exposures_corr['pheno1']
        exposures_filtered = exposures_corr.exp.to_list()
        if len(exposures_filtered) == 0:
            warnings.warn(f'No exposures found for {g2}/{p2}')
            continue
        elif len(exposures_filtered) == 1:
            warnings.warn(f'Only one exposure found for {g2}/{p2}: '+
                          f'{exposures_filtered[0]}, consider UVMR')
            continue
        
        # find instruments
        instruments = []
        exposures_corr = exposures_corr['group1'].unique()
        # SNPs identified for all exposures - cross-clumping + outcome clumped for exposure
        for i in range(len(exposures_corr)):
            g1i = exposures_corr[i]
            for j in range(i, len(exposures_corr)):
                g1j = exposures_corr[j]
                instruments.append(f'{args._in}/{g1j}_clumped_for_{g1i}_{args.pval:.0e}.txt')
            instruments.append(f'{args._in}/{g2}_clumped_for_{g1i}_{args.pval:.0e}.txt')
        
        cmd = ['Rscript mr_mvmr.r -i'] + instruments + ['--p1'] + exposures_filtered + \
            [f'--p2 {g2}/{p2} --rg {exp_corrmat} --pval {args.pval} -o {out_prefix}',force]
        submitter.add(' '.join(cmd))
    submitter.submit()
    return submitter

if __name__ == '__main__':
    from ._utils.slurm import slurm_parser
    parser = slurm_parser(description = 
      'This script batch runs Multivariable MR for groups of phenotypes')
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
    pheno_spec.add_argument('--rgp', type = float, default = -1,
                        help = 'rg p-value threshold to filter exposures, enter -1 for FDR=0.05')
    parser.add_argument('--pval', help = 'Clumping p-value threshold', default = 5e-8, type = float)
    parser.add_argument('--submodels', dest = 'sub', default = False, action = 'store_true',
                        help = 'Test all sub-models from MVMR-cML-SuSIE')
    parser.add_argument('-f','--force', dest = 'force', action = 'store_true',
                        default = False, help = 'Force overwrite')
    args = parser.parse_args()
    args.p1.sort(); args.p2.sort()
    
    import os
    for arg in ['_in','out','clump', 'gwa','rg']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    
    from ._utils import cmdhistory, path, logger
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