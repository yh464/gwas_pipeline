#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2023-07-22

This is a script for single-trait fine-mapping using polyfun-susie

Required input:
    GWAS summary statistics (scans directory for files)
'''

def main(args):
    if args.force: force = '-f'
    else: force = ''
    
    # make output directories
    import os
    if not os.path.isdir(args.out): os.mkdir(args.out)
    logdir = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/logs/'
    log = open(f'{logdir}/finemap_batch.log','w')
    
    # array submitter
    from _utils.slurm import array_submitter
    submitter = array_submitter(
        name = f'finemap_{args.pheno[0]}', n_cpu = 2,
        env = 'gentoolspy',
        timeout = 60,
        debug = False
        )
    scripts_path = os.path.realpath(__file__)
    scripts_path = os.path.dirname(scripts_path)
    
    from fnmatch import fnmatch
    for x in args.pheno:
      os.chdir(args._in)
      os.chdir(x)
      
      flist = []
      for y in os.listdir():
        if fnmatch(y,'*.fastGWA') and (not fnmatch(y,'*_X.fastGWA')) and (not fnmatch(y, '*_all_chrs*')):
          flist.append(y)
      
      for y in flist:
        skip = True
        prefix = y.replace('.fastGWA','')
        
        if not os.path.isfile(f'{args.out}/{args.pheno}/{prefix}.finemap.summary'):
          skip = False
        
        if skip and (not args.force):
          print(f'{prefix} summary already generated - skipping', file = log)
          continue
        
        print(f'{prefix} submitted to slurm', file = log)
        
        submitter.add(
          f'python finemap_by_trait.py {x} -i {y} -d {args._in} -c {args.clump} -o {args.out}'+
          f' -p {args.p:.4e} --polyfun {args.polyfun} -b {args.bfile} {force}')
    submitter.submit()
    # submitter.debug()
    
if __name__ == '__main__':
    from _utils.slurm import slurm_parser
    parser = slurm_parser(
      description = 'This programme batch runs the fine-map pipeline')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
    parser.add_argument('-i','--in', dest = '_in', help = 'Directory containing all summary stats',
      default = '../gwa/')
    parser.add_argument('-c','--clump', dest = 'clump', help = 'Directory containing all clump outputs',
      default = '../clump/')
    parser.add_argument('-o', '--out', dest = 'out', help = 'output directory',
      default = '../finemap/')
    parser.add_argument('-b', '--bfile', dest = 'bfile', help = 'bed binary',
      default = '../params/bed/')
    parser.add_argument('--polyfun', help = 'directory of POLYFUN tool',
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/polyfun/')
    parser.add_argument('-p', dest = 'p', help = 'p-value', default = 3.1076e-11, type = float)
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
      default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','clump','bfile']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    
    from _utils import path, cmdhistory, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheng/%pheno_%maf.fastGWA', __file__)
    proj.add_output(args.out+'/%pheng/%pheno_%maf.finemap.summary',__file__)
    try: main(args)
    except: cmdhistory.errlog()