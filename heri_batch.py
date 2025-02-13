#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2023-07-13
Version 2: 2024-11-14

Munges GWAS summary statistics for LDSC and then estimates heritability

Preceding workflow:
    gwa_batch.py
Requires following inputs:
    GWAS summary statistics (scans directory for all files)
'''

def main(args):
    import os
    from fnmatch import fnmatch
    
    force = '-f' if args.force else ''
    
    # array submitter
    from _utils import array_submitter
    submitter = array_submitter.array_submitter(
        name = f'heri_{args.pheno[0]}',
        timeout = 10, mode = 'long',
        debug = False
        )
    
    for x in args.pheno:
      os.chdir(args._in)
      os.chdir(x)
      if not os.path.isdir(f'{args.out}/{x}'): os.mkdir(f'{args.out}/{x}')
      
      for y in os.listdir():
        if not fnmatch(y, '*.fastGWA'): continue
        if fnmatch(y, '*X.fastGWA') or fnmatch(y, '*all_chrs.fastGWA'):
            continue                               # autosomes
        prefix = y.replace('.fastGWA','')
        if os.path.isfile(f'{args.out}/{x}/{prefix}.h2.log') and not args.force: continue
        submitter.add('python '+
            f'heri_by_trait.py -i {args._in}/{x}/{y} -o {args.out}/{x}/ --ldsc {args.ldsc} {force}')
    submitter.submit()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description = 
      'This script batch runs the LDSC heritability pipeline for local phenotypes')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
    parser.add_argument('-i','--in', dest = '_in', help = 'GWA file directory',
      default = '../gwa/')
    parser.add_argument('--ldsc', dest = 'ldsc', help = 'LDSC executable directory',
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/ldsc/') # intended to be absolute
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
      default = '../gene_corr/ldsc_sumstats/')
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
      default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','ldsc']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
        
    from _utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_var('%pheng',r'.+', 'phenotype group')
    proj.add_var('%pheno',r'.+', 'phenotype')
    proj.add_var('%maf',r'[0-9.]+', 'minor allele freq') # only allows digits and decimals
    proj.add_input(args._in+'/%pheng/%reg_%maf.fastGWA', __file__)
    proj.add_output(args.out+'/%pheng/%reg_%maf.h2.log', __file__)
    try: main(args)
    except: cmdhistory.errlog()