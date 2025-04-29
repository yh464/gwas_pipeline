#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-04-28

A flexible framework to estimate genomic SEM models

Requires following inputs: 
    MUNGED GWAS summary statistics
    FULL GWAS summary statistics (optional, only for SNP-level models)

TODO:
    - information to be passed to GenomicSEM including prevalence, effective sample size (only needed for GWAS models)
    - manual model (interactive interface)
'''

def main(args):
    # find GWAS summary stats
    from _utils.path import find_gwas
    exposures = find_gwas(args.p1, dirname=args._in, ext='sumstats', long = True) # no SE for munged summary stats
    exposures_short = find_gwas(args.p1, dirname=args._in, ext='sumstats', long = False)
    outcomes = find_gwas(args.p2, dirname=args._in, ext='sumstats', long = True)
    outcomes_short = find_gwas(args.p2, dirname=args._in, ext='sumstats', long = False)
    covariates = find_gwas(args.cov, dirname=args._in, ext='sumstats', long = True)

    from _utils.slurm import array_submitter
    submitter = array_submitter(
        name = 'gsem_'+args.p1[0]+'_'+'_'.join(args.p2), env = 'gentoolsr',
        n_cpu = 4, timeout = 60 if args.gwas else 15)
    
    # tasks string
    tasks = []
    if args.common: tasks.append('--common')
    if args.efa: tasks.append('--efa'); tasks.append(f'--efa_thr {args.efa_thr}')
    if args.mdl: tasks.append('--mdl')
    if args.gwas: tasks.append('--gwas')
    if args.force: tasks.append('--force')

    # metadata
    meta = []
    for g,_ in exposures + outcomes + covariates:
        if os.path.isfile(f'{args.full}/{g}/metadata') and not f'{args.full}/{g}/metadata' in meta: 
            meta.append(f'{args.full}/{g}/metadata')
    if len(meta) > 0: tasks += ['--meta'] + meta

    if not os.path.isdir(args.out): os.system(f'mkdir -p {args.out}')

    # if modelling exposure-outcome effects, use only correlated traits
    from gcorr_plot import crosscorr_parse
    if len(outcomes) > 0: exp_corr_out = crosscorr_parse(exposures_short, outcomes_short, logdir=args.rg)
    for g2, p2 in outcomes:
        if not os.path.isdir(f'{args.out}/{g2}'): os.system(f'mkdir -p {args.out}/{g2}')

        # screen for only correlated exposures
        if args.rgp < 0:
            exposures_corr = exp_corr_out.loc[(exp_corr_out.group2==g2) & (exp_corr_out.pheno2==p2) &\
                (exp_corr_out.q < 0.05),['group1','pheno1']]
        else:
            exposures_corr = exp_corr_out.loc[(exp_corr_out.group2==g2) & (exp_corr_out.pheno2==p2) &\
                (exp_corr_out.p < args.rgp),['group1','pheno1']]
        exposures_filtered = list(zip(exposures_corr.group1.to_list(), exposures_corr.pheno1.to_list()))

        # all correlated exposures in the same model
        if args.all_exp:
            out_prefix = f'{args.out}/{g2}/{p2}.all_'+'_'.join(args.p1)
            if os.path.isfile(f'{out_prefix}_subtraction.fastGWA') and not args.force: continue
            
            cmd = ['Rscript gsem_master.r', '-i', args._in, '-o', out_prefix, '--full', args.full, '--ref', args.ref, '--ld', args.ld,
                   '--p1'] + [f'{g}/{p}' for g, p in exposures + covariates]
            cmd += ['--p2', f'{g2}/{p2}'] + tasks
            submitter.add(' '.join(cmd))
        
        # individual correlated exposures in each model
        else:
            for g1, p1 in exposures_filtered:
                out_prefix = f'{args.out}/{g2}/{p2}.{g1}_{p1}'
                if os.path.isfile(f'{out_prefix}_subtraction.fastGWA') and not args.force: continue
                cmd = ['Rscript gsem_master.r', '-i', args._in, '-o', out_prefix, '--full', args.full, '--ref', args.ref, '--ld', args.ld,
                       '--p1', f'{g1}/{p1}']
                cmd += [f'{g}/{p}' for g, p in covariates]
                cmd += ['--p2', f'{g2}/{p2}'] + tasks
                submitter.add(' '.join(cmd))
            
    if len(outcomes) == 0:
        if args.all_exp:
            outdir = f'{args.out}/'+'_'.join(args.p1)
            if not os.path.isdir(outdir): os.system(f'mkdir -p {outdir}')
            out_prefix = f'{outdir}/{"_".join(args.p1)}_all'
            if not os.path.isfile(f'{out_prefix}.F1.fastGWA') or args.force:
                cmd = ['Rscript gsem_master.r', '-i', args._in, '-o', out_prefix, '--full', args.full, '--ref', args.ref, '--ld', args.ld,
                    '--p1'] + [f'{g}/{p}' for g, p in exposures] + tasks + ['--meta'] + meta
                submitter.add(' '.join(cmd))
        else:
            for g1, p1s in exposures_short:
                outdir = f'{args.out}/{g1}'; 
                if not os.path.isdir(outdir): os.system(f'mkdir -p {outdir}')
                cmd = ['Rscript gsem_master.r', '-i', args._in, '-o', f'{outdir}/all', '--full', args.full, '--ref', args.ref, '--ld', args.ld,
                       '--p1'] + [f'{g1}/{p}' for p in p1s] + tasks
                submitter.add(' '.join(cmd))
    
    submitter.submit()
    return submitter

if __name__ == '__main__':
    from _utils.slurm import slurm_parser
    parser = slurm_parser(description = 'A flexible framework to estimate genomic SEM models')
    path = parser.add_argument_group('Path specifications')
    path.add_argument('-i','--in', dest = '_in', help = 'MUNGED GWAS summary statistics',
        default = '../gcorr/ldsc_sumstats')
    path.add_argument('--full', help = 'FULL GWAS summary statistics (optional, only for SNP-level models)',
        default = '../gwa')
    path.add_argument('-o', '--out', help = 'Output directory', default = '../gsem')
    path.add_argument('-rg', dest = 'rg', help = 'Directory to rg log files, required for causal and subtraction',
        default = '../gcorr/rglog')
    path.add_argument('--ref', help = 'Reference file for SNP variance estimation', default = 
        '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/ldsc_for_gsem/ref.1000G.txt')
    path.add_argument('--ld', help = 'LD reference panel', default = 
        '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/ldsc/baseline')
    
    pheno = parser.add_argument_group('Phenotype specifications')
    pheno.add_argument('-p1', help = 'Exposure, scans directory', nargs = '+', default = [])
    pheno.add_argument('-p2', help = 'Outcome, scans directory', nargs = '*', default = [])
    pheno.add_argument('--cov', help = 'Covariates, scans directory', nargs = '*', default = ['demographics'])
    pheno.add_argument('--all_exp', help = 'Run common/causal/subtraction model for all exposures', 
        action = 'store_true', default = False)
    pheno.add_argument('--rgp', type = float, default = -1,
        help = 'rg p-value threshold to filter exposures, enter -1 for FDR=0.05')

    tasks = parser.add_argument_group('Analyses specifications')
    tasks.add_argument('--common', help = 'Common factor model', action = 'store_true', default = False)
    tasks.add_argument('--efa', help = 'Exploratory factor analysis', action = 'store_true', default = False)
    tasks.add_argument('--efa_thr', help = 'loading threshold to keep in a factor', default = 0.3, type = float)
    tasks.add_argument('--mdl', help = 'Causal and subtraction model', action = 'store_true', default = False)
    tasks.add_argument('--manual', help = 'Manual model', action = 'store_true', default = False)
    tasks.add_argument('--gwas', help = 'Output GWAS summary statistics', action = 'store_true', default = False)

    parser.add_argument('-f','--force', help = 'Force overwrite', action = 'store_true', default = False)
    args = parser.parse_args()

    import os
    for arg in ['_in', 'out', 'full', 'ref', 'ld', 'rg']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')

    from _utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_input(f'{args._in}/*', __file__)
    proj.add_input(f'{args.full}/*', __file__)
    proj.add_output(f'{args.out}/*',__file__)
    
    try: main(args)
    except: cmdhistory.errlog()