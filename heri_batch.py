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
  from _utils.plugins.logparser import parse_h2_log
  import numpy as np
  
  force = '-f' if args.force else ''
  
  # array submitter
  from _utils.slurm import array_submitter
  submitter = array_submitter(name = f'heri_{args.pheno[0]}', timeout = 10, env = args.ldsc, partition = 'sapphire')
  
  from _utils.path import find_gwas
  pheno = find_gwas(args.pheno, dirname = args._in, ext = 'fastGWA', long = True)

  for g, p in pheno:
    os.makedirs(f'{args.out}/{g}/', exist_ok = True)
    out_prefix = f'{args.out}/{g}/{p}'

    cmds = []
    sumstats_file = f'{args._in}/{g}/{p}.fastGWA'
    munged_file = f'{out_prefix}.sumstats'
    h2_log = f'{out_prefix}.h2.log'

    if args.force or (not os.path.isfile(munged_file)):
        hdr = open(sumstats_file).readline().strip().split()
        if 'OR' in hdr: ss = 'OR,1'
        elif 'BETA' in hdr: ss = 'BETA,0'
        elif 'Z' in hdr: ss = 'Z,0'
        else: raise ValueError('No valid summary statistics found in the input file')

        cmds.append(f'python {args.ldsc}/munge_sumstats.py --sumstats {sumstats_file} '+ \
                    f'--merge-alleles {args.ldsc}/ukb_merge_ldscore.txt '+
                    f'--signed-sumstats {ss} '+
                    f'--out {out_prefix} --chunksize 50000')
    
    if args.force or (not os.path.isfile(h2_log)):
        cmds.append(f'python {args.ldsc}/ldsc.py '+
          f'--ref-ld-chr {args.ldsc}/baseline/ --w-ld-chr {args.ldsc}/baseline/ '+
          f'--h2 {munged_file} '+
          f'--out {out_prefix}.h2')
    if len(cmds) > 0: submitter.add(*cmds)
  submitter.submit()

if __name__ == '__main__':
    from _utils.slurm import slurm_parser

    parser = slurm_parser(description = 
      'This script batch runs the LDSC heritability pipeline for local phenotypes')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
    parser.add_argument('-i','--in', dest = '_in', help = 'GWA file directory',
      default = '../gwa/')
    parser.add_argument('--ldsc', dest = 'ldsc', help = 'LDSC executable directory',
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/ldsc/') # intended to be absolute
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
      default = '../gcorr/ldsc_sumstats/')
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
      default = False, action = 'store_true')
    
    args = parser.parse_args()
    import os
    for arg in ['_in','out','ldsc']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
        
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