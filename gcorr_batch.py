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


def main(args):
    from logparser import parse_rg_log
    if not os.path.isdir(args.out): os.system(f'mkdir -p {args.out}')
    from subprocess import check_output
    
    scripts_path = os.path.dirname(os.path.realpath(__file__))
    
    # scans directories to include sumstats
    from _utils.path import find_gwas, pair_gwas
    gwa1 = find_gwas(*args.p1, dirname = args._in, ext = 'sumstats', long = False)
    gwa2 = find_gwas(*args.p2, dirname = args._in, ext = 'sumstats', long = False)
    pairwise = pair_gwas(gwa1, gwa2)
    
    # array submitter
    timeout = int(max([len(x) for _,x in (gwa1+gwa2)]+[45])/12) # each phenotype takes ~5 seconds
    from _utils.slurm import array_submitter
    submitter = array_submitter(
        name = f'gcorr_{args.p1[0]}',
        timeout = timeout, mode = 'long', wd = args._in,
        # debug = True
        )
    
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
                sumstats = [f'{g1}/{p1}.sumstats'] + \
                    [f'{g2}/{p2}.sumstats' for p2 in na_p2s]
                sumstats = ','.join(sumstats)
                submitter.add(
                    f'bash {scripts_path}/ldsc_master.sh ldsc.py '+
                    f'--ref-ld-chr {args.ldsc}/baseline/ --w-ld-chr {args.ldsc}/baseline/ '+
                    f'--rg {sumstats} --out {out_noint_rg[:-4]} --no-intercept')
        
            if os.path.isfile(out_rg) and (not args.force): continue
            sumstats = [f'{g1}/{p1}.sumstats']
            for p2 in p2s:
                if not p2 in na_p2s:
                    sumstats.append(f'{g2}/{p2}.sumstats')
            sumstats = ','.join(sumstats)
            submitter.add(
                f'bash {scripts_path}/ldsc_master.sh ldsc.py '+
                f'--ref-ld-chr {args.ldsc}/baseline/ --w-ld-chr {args.ldsc}/baseline/ '+
                f'--rg {sumstats} --out {out_rg[:-4]}')
    
    submitter.submit()
    
if __name__ == '__main__':
    import argparse
    from _utils.slurm import parser_config
    parser = argparse.ArgumentParser(description = 'This script estimates genetic cross-correlations')
    parser.add_argument('-p1', help = 'First group of phenotypes to correlate', 
                        nargs = '*', default = [])
    parser.add_argument('-p2', nargs = '*', default = [], 
        help = 'Second group of phenotypes to correlate, leave blank to calculate '+
        'pairwise correlations between all phenotypes in p1')
    parser.add_argument('-i','--in', dest = '_in', help = 'GWA file directory',
        default = '../gcorr/ldsc_sumstats/')
    parser.add_argument('--ldsc', dest = 'ldsc', help = 'LDSC executable directory',
        default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/ldsc/') # intended to be absolute
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
        default = '../gcorr/rglog/')
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
        default = False, action = 'store_true')
    parser = parser_config(parser)
    args = parser.parse_args()
    import os
    for arg in ['_in','out','ldsc']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
    from _utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheng_%pheno_%maf.sumstats', __file__)
    proj.add_output(args.out+'/%pheng_%pheno_%maf.%pheng_%pheno_%maf.rg.log', __file__)
    try: main(args)
    except: cmdhistory.errlog()