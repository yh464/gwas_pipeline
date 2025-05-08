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
    if not os.path.isdir(tmpdir): os.mkdir(tmpdir)
    
    # array submitter
    timeout = 15 if args.p < 1e-8 else 40
    from _utils.slurm import array_submitter
    submitter = array_submitter(
      name = f'clump_{args.pheno[0]}_{args.p:.0e}',
      timeout = timeout, mode = 'long')
    
    from _utils.path import find_gwas
    pheno = find_gwas(args.pheno, dirname = args._in, ext = 'fastGWA', long = True)

    # directory management
    for g,p in pheno:
      if not os.path.isdir(f'{args.out}/{g}'): os.system(f'mkdir -p {args.out}/{g}') # creates output folder
      print(f'{g}/{p}')
      out_fname = f'{args.out}/{g}/{p}_{args.p:.0e}.clumped'
      if os.path.isfile(out_fname) and (not args.force): continue
      submitter.add(
        f'python gwa_clump.py --in {args._in}/{g}/{p}.fastGWA -b {args.bfile} --plink {args.plink} '+
        f'-p {args.p} -o {args.out}/{g} {force}')
    submitter.submit()
    
if __name__ == '__main__':
    from _utils.slurm import slurm_parser
    parser = slurm_parser(description='This programme uses PLINK1.9'+
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