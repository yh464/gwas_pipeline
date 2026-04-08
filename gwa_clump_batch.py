#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2024-12-03

Clumps independent loci from a fastGWA format file

Requires following inputs: 
    GWAS summary statistics (scans directory for files)
'''

import os
from _utils import cmdhistory, path, logger
proj = path.project()
from _utils.logger import logger
log = logger()

def main(args):
    if args.force: force = '-f'
    else: force = ''
    
    # temp and log
    tmpdir = '/home/yh464/rds/rds-rb643-ukbiobank2/Data_Users/yh464/temp/'
    if not os.path.isdir(tmpdir): os.mkdir(tmpdir)
    
    # array submitter
    timeout = 15 if args.pval < 1e-8 else 40
    from _utils.slurm import array_submitter
    submitter = array_submitter(
      name = f'clump_{args.pheno[0]}_{args.pval:.0e}',
      timeout = timeout)
    
    # path specification
    proj.register('clump', 'clump/$group/$pheno_5e-8.clumped')
    for pval in args.pval:
      if pval != 5e-8: 
        proj.register(f'clump_{pval:.0e}', f'clump/$group/$pheno_{pval:.0e}.clumped')
    pheno = proj.find_gwas(args.pheno, long = True)

    # run clumping
    for g,p in pheno:
      gwa = proj.to_pathname('gwa', group = g, pheno = p)
      for pval in args.pval:
        out_prefix = proj.to_pathname('clump', group = g, pheno = p) if pval == 5e-8 else \
          proj.to_pathname(f'clump_{pval:.0e}', group = g, pheno = p)
        out_prefix = out_prefix.replace(f'_{pval:.0e}.clumped', '')
        out_fname = f'{args.out}/{g}/{p}_{pval:.0e}.clumped'
        if os.path.isfile(out_fname) and (not args.force): continue
        submitter.add(
          f'python gwa_clump.py --in {gwa} -b {args.bfile} --plink {args.plink} '+
          f'-p {pval} -o {out_prefix} {force}')
    submitter.submit()
    return submitter
    
if __name__ == '__main__':
    from _utils.slurm import slurm_parser
    parser = slurm_parser(description='This programme uses PLINK1.9'+
      ' to clump the GWAS output, identifying independent SNPs')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
    parser.add_argument('--plink', dest = 'plink', help = 'Path to PLINK *1.9* executable', 
      default = '/home/yh464/rds/rds-rb643-ukbiobank2/Data_Genetics/plink')
    parser.add_argument('-b','--bfile', dest = 'bfile', help = 'BED file list',
      default = '../params/bed')
    parser.add_argument('-p', '--pval',help = 'p-value threshold',
      default = [5e-8], type = float, nargs = '*')
    parser.add_argument('-f','--force', dest = 'force', help = 'Force output',
      default = False, action = 'store_true')
    args = parser.parse_args()
    args.bfile = os.path.realpath(args.bfile)

    from _utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    try: main(args)
    except: cmdhistory.errlog()