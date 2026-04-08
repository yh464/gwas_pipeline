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

import os
from _utils import cmdhistory, path, logger
proj = path.project()

def main(args):
  # array submitter
  from _utils.slurm import array_submitter
  submitter = array_submitter(name = f'heri_{args.pheno[0]}', timeout = 10 if not args.complete else 60, env = args.ldsc, partition = 'sapphire')
  
  # path specification
  proj.register('ldsc_sumstats', 'gcorr/ldsc_sumstats/$group/$pheno.sumstats')
  if args.complete:
    proj.register('ldsc_sumstats_complete', 'gcorr/ldsc_sumstats/$group_complete/$pheno.sumstats')
  pheno = proj.find_gwas(args.pheno, long = True)

  for g, p in pheno:
    out_prefix = proj.to_pathname('ldsc_sumstats', group = g, pheno = p) if not args.complete else \
      proj.to_pathname('ldsc_sumstats_complete', group = g, pheno = p)
    out_prefix = out_prefix.rsplit('.sumstats', 1)[0] # remove .sumstats suffix as LDSC will add it back
    sumstats_file = proj.to_pathname('gwa', group = g, pheno = p)
    munged_file = f'{out_prefix}.sumstats'
    h2_log = f'{out_prefix}.h2.log'
    cmds = []

    if args.force or (not os.path.isfile(munged_file)):
      hdr = open(sumstats_file).readline().strip().split()
      if 'OR' in hdr: ss = 'OR,1'
      elif 'BETA' in hdr: ss = 'BETA,0'
      elif 'Z' in hdr: ss = 'Z,0'
      else: raise ValueError('No valid summary statistics found in the input file')

      cmds.append(f'python {args.ldsc}/munge_sumstats.py --sumstats {sumstats_file} '+ \
                  (f'--merge-alleles {args.ldsc}/ukb_snp_info.txt ' if args.complete else \
                  f'--merge-alleles {args.ldsc}/ukb_merge_ldscore.txt ')+
                  f'--signed-sumstats {ss} '+
                  f'--out {out_prefix} --chunksize 50000')
      cmds.append(f'if [ -f {out_prefix}.sumstats.gz ]; then gunzip -f {out_prefix}.sumstats.gz; fi')
    
    if args.force or (not os.path.isfile(h2_log)):
      cmds.append(f'python {args.ldsc}/ldsc.py '+
        f'--ref-ld-chr {args.ldsc}/baseline/ --w-ld-chr {args.ldsc}/baseline/ '+
        f'--h2 {munged_file} '+
        f'--out {out_prefix}.h2')
    if len(cmds) > 0: submitter.add(*cmds)
  submitter.submit()
  return submitter

def api(**kwargs):
  from _utils.gadgets import namespace
  args = namespace(**kwargs)
  return main(args)

if __name__ == '__main__':
    from _utils.slurm import slurm_parser
    parser = slurm_parser(description = 
      'This script batch runs the LDSC heritability pipeline for local phenotypes')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
    parser.add_argument('--ldsc', dest = 'ldsc', help = 'LDSC executable directory',
      default = '/home/yh464/rds/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/ldsc/') # intended to be absolute
    # remove the input and output options as these files are expected to be in standard locations
    parser.add_argument('-c','--complete', action = 'store_true', help = 'Merge with the complete set of SNPs in UKB, instead of HapMap3')
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
      default = False, action = 'store_true')
    args = parser.parse_args()
    args.ldsc = os.path.realpath(args.ldsc)
        
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()