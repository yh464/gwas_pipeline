#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2024-11-05
Version 2: 2025-01-21

Batch submits jobs for Mendelian Randomisation for all GWAS files in a directory
(usually the same group of phenotypes). Scans the entire directory for GWAS summary
stats of the same data extension.
'''
def qc(file):
    import pandas as pd
    try: 
        df = pd.read_table(file).dropna()
        return df.size >= 5
    except:
        return False

def main(args):
    import os
    from ._plugins.logparser import parse_h2_log
    from ._utils.path import find_clump, find_gwas
    from ._plugins.logparser import crosscorr_parse
    
    # array submitter
    from mr_extract_snp_batch import api
    print('Try running following command to force re-extraction of instruments')
    print(f'python mr_extract_snp_batch.py -p1 {" ".join(args.p1)} -p2 {" ".join(args.p2)} -b -i {args.gwa} -o {args.inst} -c {args.clump} -f')
    snp_submitter = api(p1 = args.p1, p2 = args.p2, bid = True, _in = args.gwa, out = args.inst, clump = args.clump)
    from ._utils.slurm import array_submitter
    submitter_main = array_submitter(name = 'mr_'+'_'.join(args.p2), env = 'gentoolsr',
        n_cpu = 3 if args.apss else 2, timeout = 7, dependency = snp_submitter)
    submitter_lcv = array_submitter(name = 'mr_lcv_'+'_'.join(args.p2), env = 'gentoolsr', n_cpu = 2, timeout = 30)
    submitter_cause = array_submitter(name = 'mr_cause_'+'_'.join(args.p2), env = 'gentoolsr',n_cpu = 3, timeout = 30)
    
    # output directory
    if not os.path.isdir(args.out): os.mkdir(args.out)

    # genetic correlations
    exposures = find_gwas(args.p1, dirname=args.gwa, ext=args.ext1, se = True)
    outcomes = find_gwas(args.p2, dirname=args.gwa, ext=args.ext2, se = True)
    exp_corr_out = crosscorr_parse(exposures, outcomes, logdir=args.rg, full = True)

    # general command args for mr_master
    force = '-f' if args.force else ''
    cmdargs = []
    if args.apss: cmdargs.append('--apss')
    meta = []
    for g, _ in exposures + outcomes:
        if os.path.isfile(f'{args.gwa}/{g}/metadata') and not f'{args.gwa}/{g}/metadata' in meta: 
            meta.append(f'{args.gwa}/{g}/metadata')
    if len(meta) > 0: cmdargs.append(f'--meta {":".join(meta)}')

    for g2, p2s in outcomes:
      # make output directory wrt g2
      if not os.path.isdir(f'{args.out}/{g2}'): os.mkdir(f'{args.out}/{g2}')

      for g1, p1s in exposures:
        for p1 in p1s:
          # find h2 log for trait 1
          h2log = f'{args.h2}/{g1}/{p1}.h2.log'
          h21, h2se1 = parse_h2_log(h2log)
          
          for p2 in p2s:
            if not os.path.isdir(f'{args.out}/{g2}/{p2}'): 
                os.mkdir(f'{args.out}/{g2}/{p2}')  
            
            # filter for rg
            rginfo = exp_corr_out.loc[(exp_corr_out.group1==g1) & (exp_corr_out.pheno1==p1) & \
                (exp_corr_out.group2==g2) & (exp_corr_out.pheno2==p2),:]
            if rginfo.shape[0] == 0: print(f'No rg for {g1}/{p1} and {g2}/{p2}'); continue
            if (args.rgp > 0 and rginfo.p.values[0] > args.rgp) or (args.rgp < 0 and rginfo.q.values[0] > 0.05):
                print(f'Correlation between {g1}/{p1} and {g2}/{p2} is not significant, skipping')
                continue
            
            # find h2 log for trait 1
            h2log = f'{args.h2}/{g2}/{p2}.h2.log'
            h22, h2se2 = parse_h2_log(h2log)
            
            # input file name specification
            try:
                gwa1 = f'{args.gwa}/{g1}/{p1}.{args.ext1}'
                clump1, pval1 = find_clump(g1, p1, args.clump, args.pval)
                clump001, _ = find_clump(g1, p1, args.clump, 0.001)
            except: print(f'{g1} missing clumped GWAS sumstats'); continue
            try: 
                gwa2 = f'{args.gwa}/{g2}/{p2}.{args.ext2}'
                clump2, pval2 = find_clump(g2, p2, args.clump, args.pval)
                clump002, _ = find_clump(g2, p2, args.clump, 0.001)
            except: print(f'{g2} missing clumped GWAS sumstats'); continue
            pval_thr = max([pval1, pval2])
            
            # find instruments
            instruments = []
            for i,_ in exposures:
                for j,_ in outcomes:
                    instruments.append(f'{args.inst}/{i}_clumped_for_{j}_{pval_thr:.0e}.txt')
                    instruments.append(f'{args.inst}/{j}_clumped_for_{i}_{pval_thr:.0e}.txt')
                    instruments.append(f'{args.inst}/{j}_clumped_for_{j}_{pval_thr:.0e}.txt')
                instruments.append(f'{args.inst}/{i}_clumped_for_{i}_{pval_thr:.0e}.txt')

            # check progress and QC output
            out_prefix = f'{args.out}/{g2}/{p2}/{g1}_{p1}_{p2}'
            fwd = f'{out_prefix}_mr_forward_results.txt'
            rev = f'{out_prefix}_mr_reverse_results.txt'
            fwd_presso = f'{out_prefix}_mr_forward_presso_results.txt'
            rev_presso = f'{out_prefix}_mr_reverse_presso_results.txt'
            
            for file in [fwd, rev, fwd_presso, rev_presso]:
                if not os.path.isfile(file): continue
                if not qc(file):
                    try: os.remove(file)
                    except: pass
                    print(f'Empty file: {file}')
            
            if not os.path.isfile(fwd) or \
                not os.path.isfile(rev) or \
                not os.path.isfile(fwd_presso) or \
                not os.path.isfile(rev_presso) or args.force:
                cmd = ['Rscript mr_master.r','--p1', f'{g1}/{p1}', '--p2', f'{g2}/{p2}',
                    '-i', ':'.join(instruments), '--c1', clump1, '--c2', clump2,
                    '--g1', gwa1, '--g2', gwa2, '--pval', str(pval_thr), '--h21', str(h21),
                    '--h2se1', str(h2se1), '--h22', str(h22), '--h2se2', str(h2se2),
                    '--gcovint', str(rginfo.gcov_int.values[0]), '--gcintse', str(rginfo.gcov_int_se.values[0]),
                    force, '--ldsc', args.ldsc,'-o', f'{args.out}/{g2}/{p2}'] + cmdargs
                submitter_main.add(' '.join(cmd))
            
            if not os.path.isfile(f'{out_prefix}_mr_lcv_results.txt') or args.force: # bidirectional
                submitter_lcv.add(
                    f'Rscript mr_lcv.r --g1 {gwa1} --g2 {gwa2} '+
                    f'-o {args.out}/{g2}/{p2} --ldsc {args.ldsc} {force}')
            
            if not os.path.isfile(f'{out_prefix}_mr_forward_cause_results.txt') or \
                not os.path.isfile(f'{out_prefix}_mr_reverse_cause_results.txt') or args.force:
                submitter_cause.add(
                    f'Rscript mr_cause.r --g1 {gwa1} --c1 {clump001} '+
                    f'--g2 {gwa2} --c2 {clump002} -o {args.out}/{g2}/{p2} {force}')
    
    submitter_main.submit()
    submitter_lcv.submit()
    submitter_cause.submit()
    
if __name__ == '__main__':
    from ._utils.slurm import slurm_parser
    parser = slurm_parser(description = 
      'This script batch runs MR for groups of phenotypes')
    path_spec = parser.add_argument_group('Path specifications')
    path_spec.add_argument('-g','--gwa', dest = 'gwa', default = '../gwa',
        help = 'input GWA directory, assumes both groups of pheno to be in the same dir')
    path_spec.add_argument('-i','--inst', dest = 'inst', help = 'Instruments from mr_extract_snp',
        default = '../mr/instruments')
    path_spec.add_argument('-c','--clump', dest = 'clump', help = 'Directory of clumping files',
        default = '../clump')
    path_spec.add_argument('-h2', dest = 'h2', help = 'Directory to h2 log files',
        default = '../gcorr/ldsc_sumstats')
    path_spec.add_argument('-rg', dest = 'rg', help = 'Directory to rg log files',
        default = '../gcorr/rglog')
    path_spec.add_argument('-o','--out', dest = 'out', help = 'Output directory',
        default = '../mr')
    path_spec.add_argument('--ldsc', help = 'LD scores, for LCV regression', # intentionally absolute
        default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/ldsc_for_gsem/uk10k.l2.ldscore')
    
    pheno_spec = parser.add_argument_group('Phenotype specifications')
    pheno_spec.add_argument('-p1','--pheno1', dest = 'p1', 
        help = 'Phenotypes group 1', nargs = '*', default=[])
    pheno_spec.add_argument('-e1','--ext1', dest = 'ext1', help = 'Extension for phenotype group 1',
        default = 'fastGWA')
    pheno_spec.add_argument('-p2','--pheno2', dest = 'p2', 
        help = 'Phenotypes group 2', nargs = '*', default = [])
    pheno_spec.add_argument('-e2','--ext2', dest = 'ext2', help = 'Extension for phenotype group 2',
        default = 'fastGWA')
    pheno_spec.add_argument('--rgp', type = float, default = -1,
        help = 'rg p-value threshold to filter exposures, enter -1 for FDR=0.05')
    
    corr = parser.add_argument_group('Corrections and adjustments')
    corr.add_argument('--apss', help = 'conduct MR-APSS correction', action = 'store_true', default = False)
    
    parser.add_argument('--pval', help = 'Clumping p-value threshold', default = 5e-8, type = float)
    parser.add_argument('-f','--force', dest = 'force', action = 'store_true',
                        default = False, help = 'Force overwrite')
    args = parser.parse_args()
    
    import os
    for arg in ['gwa','out','inst','clump', 'h2','rg']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    
    from ._utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_input(f'{args.gwa}/{args.p1}/*.{args.ext1}', __file__)
    proj.add_input(f'{args.gwa}/{args.p2}/*.{args.ext2}', __file__)
    proj.add_input(f'{args.clump}/{args.p1}/*.clumped',__file__)
    proj.add_input(f'{args.clump}/{args.p2}/*.clumped',__file__)
    proj.add_output(f'{args.out}/{args.p2}/*',__file__)
    
    try: main(args)
    except: cmdhistory.errlog()