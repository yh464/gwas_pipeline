#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2024-11-29

Runs GWAS meta-analysis across datasets

Requires following inputs:
    GWAS summary statistics (put files from same dataset in same folder, scans dir)
    requires file names to be identical across datasets
    requires fastGWA format (SNP, A1, A2, AF1, BETA, P)
'''

def main(args):
    import os
    from fnmatch import fnmatch
    
    # array submitter
    from _utils.slurm import array_submitter
    submitter = array_submitter(
        name = 'gwa_meta_'+'_'.join(args.dsets),
        timeout = 10, mode = 'long',
        debug = False
        )
    
    force = ' -f' if args.force else ''
    
    os.chdir(args._in)
    # scans dirs for fastGWA files
    flist = []
    for x in args.dsets:
        for y in os.listdir(x):
            if fnmatch(y, '*.fastGWA') and not fnmatch(y, '*all_chrs*')\
                and not fnmatch(y, '*_X.fastGWA') and not y in flist:
                flist.append(y)
    
    # tmp and output dir
    if not os.path.isdir(f'{args._in}/{args.out}'): os.mkdir(f'{args._in}/{args.out}')
    
    for y in flist:
        out = f'{args._in}/{args.out}/{y}'
        cmd = 'python gwa_meta.py -i '
        for x in args.dsets:
            if os.path.isfile(f'{args._in}/{x}/{y}'):
                cmd += f'{args._in}/{x}/{y} '
        cmd += f'--metal {args.metal} --plink {args.plink} -o {out} {force}'
        submitter.add(cmd)
    submitter.submit()
    
if __name__ == '__main__':
    from _utils.slurm import slurm_parser
    parser = slurm_parser(description = 
      'This programme creates genetic correlation matrices for global phenotypes')
    parser.add_argument('dsets', nargs = '*',
        help = 'Datasets to meta-analyse, scans directories for summary stats (fastGWA format)')
    parser.add_argument('-i','--in', dest = '_in', help = 'GWA file directory',
      default = '../gwa/')
    parser.add_argument('--metal', help = 'METAL executable',
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/metal') # intended to be absolute
    parser.add_argument('--plink', help = 'PLINK 1.9 executable',
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/plink') # intended to be absolute
    parser.add_argument('-o','--out', dest = 'out', 
      help = 'output directory, relative to the --in dir')
    parser.add_argument('-f','--force',dest = 'force', help = 'force overwrite',
      default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['_in','metal']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
    from _utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    for x in args.dsets: proj.add_output(args._in+'/'+x, __file__)
    proj.add_output(args._in+'/'+args.out, __file__)
    try: main(args)
    except: cmdhistory.errlog()