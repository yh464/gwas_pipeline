#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-07-31

Utility script to extract given locus from GWAS files

Requires following inputs: 
    GWAS summary statistics (single file)
'''
def main(args = None, **kwargs):
  import os
  from _utils.gadgets import namespace
  from hashlib import sha256
  if args == None: args = namespace(**kwargs) # used for API
  
  # find clumped loci
  from _utils.path import find_gwas
  from _plugins.logparser import parse_clump
  pheno = find_gwas(args.pheno, dirname = args._in, ext = 'fastGWA', clump = True, long = True)
  _, loci = parse_clump(pheno, clump_dir = args.clump, pval = args.pval)
  loci = loci.loc[loci.P < args.pval, ['CHR', 'START', 'STOP']]
  loci['START'] -= 5e5; loci['STOP'] += 5e5
  tmpdir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/temp/coloc'
  if not os.path.isdir(tmpdir): os.system(f'mkdir -p {tmpdir}')
  tempfile = f'{tmpdir}/{sha256(repr(pheno).encode()).hexdigest()}.txt'
  loci.to_csv(tempfile, sep = '\t')

  from _utils.slurm import array_submitter
  submitter = array_submitter(name = f'finemap_extract_locus_{args.pheno[0]}', n_cpu = 2, timeout = loci.shape[0]/3)
  for g, p in pheno:
    out_file = f'{args.out}/{g}/{p}_chr{loci.CHR.iloc[-1]}_{loci.START.iloc[-1]:.0f}_{loci.STOP.iloc[-1]:.0f}.txt'
    if not os.path.isfile(out_file) or args.force:
      cmd = ['python', 'gwa_extract_locus.py', f'{g}/{p}', 
             '-i', args._in, '--loci', tempfile, '-o', f'{args.out}/loci']
      submitter.add(' '.join(cmd))
  submitter.submit()
  return submitter


if __name__ == '__main__':
  from _utils.slurm import slurm_parser
  parser = slurm_parser(description = 'Extracts summary statistics for significant loci for colocalisation analysis')
  parser.add_argument('pheno', nargs = '*', default = [], help = 'Phenotypes to colocalise')
  parser.add_argument('-i','--in', dest = '_in', help = 'Input directory of GWAS summary stats',
    default = '../gwa')
  parser.add_argument('-c','--clump', help = 'Directory of clumping output',
    default = '../clump')
  parser.add_argument('-p','--pval', type = float, default = 5e-8, help = 'P-value threshold for clumping')
  parser.add_argument('-o','--out', help = 'Output directory',
    default = '../coloc/')
  parser.add_argument('-f','--force', action = 'store_true', default = False, help = 'Force overwrite')
  args = parser.parse_args()

  import os
  for arg in ['_in','out','clump']:
    setattr(args, arg, os.path.realpath(getattr(args, arg)))
  
  from _utils import cmdhistory, logger
  logger.splash(args)
  cmdhistory.log()
  try: main(args)
  except: cmdhistory.errlog()