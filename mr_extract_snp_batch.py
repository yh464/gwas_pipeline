#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-03-24

A tool to extract instruments for MR analysis before mr_batch and mr_mvmr_batch
'''
    
def main(args):
    import os
    import pandas as pd
    from _utils.path import find_gwas, find_clump
    
    # array submitter
    from _utils.slurm import array_submitter
    submitter = array_submitter(name = f'mr_instruments_'+'.'.join((args.p1+args.p2)[:2]),
        n_cpu = 1, timeout = 60)
    
    # output directory
    if not os.path.isdir(args.out): os.system(f'mkdir -p {args.out}')
    tmpdir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/temp/clump_snp_list'
    if not os.path.isdir(tmpdir): os.system(f'mkdir -p {tmpdir}')
    
    force = ' -f' if args.force else ''
    
    # parse list of outcomes and exposures
    if not args.bid: exposures = find_gwas(args.p1)
    else: exposures = find_gwas(args.p1 + args.p2)
    all_pheno = find_gwas(args.p1 + args.p2)
    
    print('Trying to find genetic instruments for all exposure phenotypes')
    for expg, expp in exposures:
        print(f'    {expg}: {len(expp)} phenotypes')
        all_snps = []
        # read trait-wise clump files
        for x in expp:
            clump,_ = find_clump(expg, x, args.clump, args.pval)
            all_snps.append(pd.read_table(clump, sep = '\\s+', usecols = ['SNP']))
        all_snps = pd.concat(all_snps)['SNP'].unique()
        temp_snps = f'{tmpdir}/{expg}_{args.pval:.0e}.txt'
        with open(temp_snps,'w') as file:
            for x in all_snps: print(x, file=file)
            file.close()
        
        for p,_ in all_pheno:
            out_file = f'{args.out}/{p}_clumped_for_{expg}_{args.pval:.0e}.txt'
            if os.path.isfile(out_file) and not args.force: continue
            cmd = f'python gwa_extract_snp.py {temp_snps} -p {p} -o {out_file} {force}'
            submitter.add(cmd)
    submitter.submit()
    
    return submitter

def api(**kwargs):
    import os
    class c(): pass
    args = c()
    # defaults
    for key, value in zip(['_in','clump','out','p1','p2','pval','bid','force'],
        ['../gwa','../clump','../mr/instruments', [],[],5e-8,False, False]):
        setattr(args, key, value)
    for key, value in kwargs.items(): setattr(args, key, value)
    for arg in ['_in','out','clump']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    submitter = main(args)
    return submitter

if __name__ == '__main__':
    from _utils.slurm import slurm_parser
    parser = slurm_parser(description = 
      'This script batch extracts instruments for MR for groups of phenotypes')
    path_spec = parser.add_argument_group('Path specifications')
    path_spec.add_argument('-i','--in', dest = '_in', 
                        help = 'input GWA directory, assumes both groups of pheno to be in the same dir',
                        default = '../gwa')
    path_spec.add_argument('-c','--clump', dest = 'clump', help = 'Directory of clumping files',
                        default = '../clump')
    path_spec.add_argument('-o','--out', dest = 'out', help = 'Output directory',
                        default = '../mr/instruments')
    pheno_spec = parser.add_argument_group('Phenotype specifications')
    pheno_spec.add_argument('-p1','--pheno1', dest = 'p1', 
                        help = 'Exposure', nargs = '*', default = [])
    pheno_spec.add_argument('-p2','--pheno2', dest = 'p2', 
                        help = 'Outcome', nargs = '*', default = [])
    parser.add_argument('--pval', help = 'Clumping p-value threshold', default = 5e-8, type = float)
    parser.add_argument('-b','--bidirectional', default = False, action = 'store_true',
                        dest = 'bid', help = 'Also extract for reverse direction')
    parser.add_argument('-f','--force', dest = 'force', action = 'store_true',
                        default = False, help = 'Force overwrite')
    args = parser.parse_args()
    if len(args.p2) == 0: args.p2 = args.p1; args.bid = False
    args.p1.sort(); args.p2.sort()
    
    import os
    for arg in ['_in','out','clump']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
    from _utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_input(f'{args._in}/{args.p1}/*.fastGWA', __file__)
    proj.add_input(f'{args._in}/{args.p2}/*.fastGWA', __file__)
    proj.add_input(f'{args.clump}/{args.p1}/*.clumped',__file__)
    proj.add_input(f'{args.clump}/{args.p2}/*.clumped',__file__)
    proj.add_output(f'{args.out}/*',__file__)
    
    try: main(args)
    except: cmdhistory.errlog()