#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2024-12-03

Clumps independent loci from a fastGWA format file

Requires following inputs: 
    GWAS summary statistics (scans directory for files)
'''

def main(args):
    import os
    from fnmatch import fnmatch
    
    if args.force: force = '-f'
    else: force = ''
    
    # temp and log
    tmpdir = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/temp/'
    logdir = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/logs/'
    if not os.path.isdir(tmpdir): os.mkdir(tmpdir)
    if not os.path.isdir(logdir): os.mkdir(logdir)
    log = open(f'{logdir}/gwa_clump.log','w')                                       # log file
    
    # array submitter
    timeout = 5 if args.p < 1e-8 else 20
    from _utils import array_submitter
    submitter = array_submitter.array_submitter(
        name = f'clump_{args.pheno[0]}_{args.p:.0e}',
        timeout = timeout,
        debug = False
        )
    
    # directory management
    for x in args.pheno:
      if not os.path.isdir(f'{args.out}/{x}'): os.system(f'mkdir -p {args.out}/{x}') # creates output folder
      os.chdir(args._in)
      os.chdir(x)
      
      # filter out required GWA files
      flist = []
      for y in os.listdir():
          if fnmatch(y, '*_X.fastGWA'): continue
        # if fnmatch(y, '*_all_chrs.fastGWA') or fnmatch(y, '*gwama.fastGWA'):
          if fnmatch(y, '*.fastGWA') or fnmatch(y, '*.txt'): flist.append(y)
      print(flist, file = log)
      
      os.chdir(args.out)
      os.chdir(x)
      for y in flist:
        out_fname = y.replace('fastGWA','clumped')
        if os.path.isfile(out_fname) and (not args.force): continue
        scripts_path = os.path.realpath(__file__)
        scripts_path = os.path.dirname(scripts_path)
        submitter.add(
          # f'bash {scripts_path}/pymaster.sh '+
          f'python gwa_clump.py --in {args._in}/{x}/ --file {y} -b {args.bfile} --plink {args.plink} '+
          f'-p {args.p} -o {args.out}/{x}/ {force}')
    submitter.submit()
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='This programme uses PLINK1.9'+
      ' to clump the GWAS output, identifying independent SNPs')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*',
      default=['deg_local','degi_local','degc_local','clu_local','eff_local','mpl_local'])
    parser.add_argument('-i','--in', dest = '_in', help = 'Input directory',
      default = '../gwa/')
    parser.add_argument('--plink', dest = 'plink', help = 'Path to PLINK *1.9* executable', 
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Genetics/plink')
    parser.add_argument('-b','--bfile', dest = 'bfile', help = 'BED file list',
      default = '../params/bed_files_ukb.txt')
    parser.add_argument('-o','--out', dest = 'out', help = 'Output directory',
      default = '../clump/')
    parser.add_argument('-p',help = 'p-value threshold',
      default = 5e-8, type = float) # or 3.1076e-11, or 1e-6
    parser.add_argument('-f','--force', dest = 'force', help = 'Force output',
      default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','bfile']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')

    from _utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_var('%pheng',r'.+', 'phenotype group')
    proj.add_var('%pheno',r'.+', 'phenotype')
    proj.add_var('%maf',r'[0-9.]+', 'minor allele freq') # only allows digits and decimals
    proj.add_input(args._in+'/%pheng/%pheno_%maf.fastGWA', __file__)
    proj.add_input(args.out+'/%pheng/%pheno_%maf.clumped', __file__)
    try: main(args)
    except: cmdhistory.errlog()