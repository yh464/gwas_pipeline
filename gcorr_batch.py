#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2024-11-14
Version 2: 2025-04-09

A simplified script to conduct genetic correlation between groups of phenotypes

Requires following inputs: 
    GWAS summary statistics (scans directory for all files)
Outputs:
    rg log between one phenotype and all phenotypes of a group
'''

import os
from _utils import cmdhistory, path, logger
proj = path.project()

def main(args):
    from _utils.plugins.logparser import parse_rg_log
    from subprocess import check_output
    
    # check progress in ldsc formatting
    from _utils.path import pair_gwas
    gwa1 = proj.find_gwas(*args.p1, long = False)
    gwa2 = proj.find_gwas(*args.p2, long = False)
    pairwise = pair_gwas(gwa1, gwa2)
    sumstats_ftype = 'ldsc_sumstats_complete' if args.complete else 'ldsc_sumstats'
    to_munge = []
    for g, ps in gwa1 + gwa2:
        for p in ps:
            sumstats = proj.to_pathname(sumstats_ftype, group = g, pheno = p)
            if not os.path.isfile(sumstats): to_munge.append(f'{g}/{p}')
    if len(to_munge) > 0:
        from heri_batch import api
        dep = api(pheno = to_munge, ldsc = args.ldsc, complete = args.complete)
    else: dep = []
    
    # input and output directory
    if not os.path.isdir(args.out): os.system(f'mkdir -p {args.out}')
    wd = proj.project_root + '/' + proj.config[sumstats_ftype].split('/$group')[0]

    # array submitter
    timeout = int(max([len(x) for _,x in (gwa1+gwa2)]+[45])/12) # each phenotype takes ~5 seconds
    from _utils.slurm import array_submitter
    submitter = array_submitter(name = f'gcorr_{gwa1[0][0]}',timeout = timeout, wd = wd, env = args.ldsc, dependency = dep)
    
    for g1, p1s, g2, p2s in pairwise:
        # p1s means list of <pheno1>s in group1
        # output data structure: {args.out}/<group1>.<group2>/<group1>_<pheno1>.<group2>.rg.log
        # where group1 <= group2
        if g1 > g2: g1, g2, p1s, p2s = g2, g1, p2s, p1s
        if not os.path.isdir(f'{args.out}/{g1}.{g2}'): os.mkdir(f'{args.out}/{g1}.{g2}')
        
        for p1 in p1s:
            if g1 == g2:
                p2s = p1s[p1s.index(p1):] # not p1s.index(p1) + 1 lest it throws an error on last element
                p2s.remove(p1)
            if len(p2s) == 0: continue
            out_rg = f'{args.out}/{g1}.{g2}/{g1}_{p1}.{g2}.rg.log'
            
            # QC out_rg file to identify NA correlations
            na_p2s = []
            if os.path.isfile(out_rg):
                # check that analysis has finished
                eof = check_output(['tail',out_rg,'-n','2']).decode()
                line = eof.split('\n')[0]
                if line.find('Analysis finished') == -1: os.remove(out_rg)
            if os.path.isfile(out_rg):
                # check for NA correlations
                all_rg = parse_rg_log(out_rg)
                if all_rg.shape[0] == 0: os.remove(out_rg)
                # check for updates to the sumstats
                if (len(all_rg.pheno2.unique()) < len(p2s)): os.remove(out_rg)
                na_p2s = all_rg.loc[all_rg.rg.isna(),'pheno2'].tolist()
                del all_rg
            
            # for NA correlations, run with constrained intercepts
            out_noint_rg = out_rg.replace('.rg.log','.noint.rg.log')
            if len(na_p2s) > 0 and (not os.path.isfile(out_noint_rg) or args.force):
                sumstats = [proj.to_pathname(sumstats_ftype, group = g1, pheno = p1).replace(f'{wd}/','')] + \
                    [proj.to_pathname(sumstats_ftype, group = g2, pheno = p2).replace(f'{wd}/','') for p2 in na_p2s]
                sumstats = ','.join(sumstats)
                submitter.add(
                    f'python {args.ldsc}/ldsc.py '+
                    f'--ref-ld-chr {args.ldsc}/baseline/ --w-ld-chr {args.ldsc}/baseline/ '+
                    f'--rg {sumstats} --out {out_noint_rg[:-4]} --no-intercept')
        
            if os.path.isfile(out_rg) and (not args.force): continue
            sumstats = [proj.to_pathname(sumstats_ftype, group = g1, pheno = p1).replace(f'{wd}/','')] + [
                proj.to_pathname(sumstats_ftype, group = g2, pheno = p2).replace(f'{wd}/','') for p2 in p2s if p2 not in na_p2s]
            if len(sumstats) < 2: continue
            sumstats = ','.join(sumstats)
            submitter.add(
                f'python {args.ldsc}/ldsc.py '+
                f'--ref-ld-chr {args.ldsc}/baseline/ --w-ld-chr {args.ldsc}/baseline/ '+
                f'--rg {sumstats} --out {out_rg[:-4]}')
    
    submitter.submit()
    
if __name__ == '__main__':
    from _utils.slurm import slurm_parser
    parser = slurm_parser(description = 'This script estimates genetic cross-correlations')
    parser.add_argument('-p1', help = 'First group of phenotypes to correlate', 
                        nargs = '*', default = [])
    parser.add_argument('-p2', nargs = '*', default = [], 
        help = 'Second group of phenotypes to correlate, leave blank to calculate '+
        'pairwise correlations between all phenotypes in p1')
    # removed input as it should be specified in the project path spec
    parser.add_argument('-c','--complete', action = 'store_true', 
        help = 'Merge with the complete set of SNPs in UKB, instead of HapMap3')
    parser.add_argument('--ldsc', dest = 'ldsc', help = 'LDSC executable directory',
        default = '/home/yh464/rds/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/ldsc/') # intended to be absolute
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
        default = '../gcorr/rglog/')
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
        default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['out','ldsc']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    try: main(args)
    except: cmdhistory.errlog()