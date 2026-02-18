#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2023-07-04
Version 2: 2025-02-24

A script to batch run GWAS by fastGWA-mlm

Requires following inputs: 
    phenotype file in FID IID *** format
    PLINK bed binaries
    covariates files in FID IID *** format
'''
from _utils.logger import logger
log = logger()
import os, fnmatch
import pandas as pd

def main(args):
  # array submitter
  from _utils.slurm import array_submitter
  submitter = array_submitter(name = 'gwa_'+ '_'.join(args.pheno),timeout = 90)
  
  # general args
  force = '-f' if args.force else ''
  xchr = '' if args.xchr else '--nox'
  if len(args.extract) > 0:
    if os.path.isfile(args.extract[0]): extract = f'--extract {os.path.realpath(args.extract[0])} '
    else:
      snp_file = f'{submitter.tmpdir}/snps_to_extract.txt'
      with open(snp_file, 'w') as f:
        for snp in args.extract:
          print(snp, file = f)
      extract = f'--extract {snp_file} '
  else: extract = ''

  # locate phenotype file
  for pheno in args.pheno:
    flist = []
    for f in os.listdir(args._in):
      if fnmatch.fnmatch(f,f'*{pheno}*.txt') and not(os.path.isdir(f)):       # search for all files matching pheno
        flist.append(f'{args._in}/{f}')
    if len(flist) != 1: log.warn(f'Please give only ONE phenotype file for {pheno}'); continue
    f = flist[0]

    # check validity of the phenotype file
    hdr = open(f).readline().replace('\n','').split()
    if hdr[0] != 'FID' or hdr[1] != 'IID':
      log.warn('Phenotype file should be in the format: FID IID *pheno')
      continue
    
    # create output folder
    outdir = f'{args.out}/{os.path.basename(f)}/'.replace('.txt','')
    log.log(outdir)
    if not os.path.isdir(outdir):
      os.system(f'mkdir -p {outdir}')                                              # this also generates args.out
    
    # phenotypes to be analysed
    c = hdr[2:]
    log.log(f'Following traits are to be GWA-analysed for {pheno}:')
    for i in c: log.log(f'    {i}')
    
    # for each phenotype
    for i, trait in enumerate(c):
      mpheno = i+1
      out_fname = outdir + trait
      # check existing files
      if os.path.isfile(f'{out_fname}.fastGWA') and not args.force:
        log.log(f'Trait already analysed for: {trait}')
        continue
      
      submitter.add(
        f'python gwa_by_trait.py -i {f} -o {out_fname} --mpheno {mpheno} --dcov {args.dcov} '+
        f'--qcov {args.qcov} --bed {args.bed} --grm {args.grm} --gcta {args.gcta} --maf {args.maf} '+
        f'--keep {args.keep} {xchr} --xbed {args.xbed} {extract} {force}'
        )
  submitter.submit()

if __name__ == '__main__':
  from _utils.slurm import slurm_parser
  parser = slurm_parser(description=
    'This programme runs GWA for any phenotype given as the 1st positional argument')
  parser.add_argument('pheno', nargs = '+', help = 'Phenotype file in TXT format - please supply ONLY ONE')
  
  io = parser.add_argument_group(title = 'input and output options')
  io.add_argument('-i','--in', dest = '_in', help = 'Phenotype directory',
    default = '../pheno/ukb/')
  io.add_argument('-o','--out',dest  = 'out', help = 'Output directory',
    default = '../gwa/')
  io.add_argument('--dcov',dest = 'dcov', help = 'DISCRETE covariance file',
    default = '../params/ukb_dcov.txt')
  io.add_argument('--qcov',dest = 'qcov', help = 'QUANTITATIVE covariance file',
    default = '../params/ukb_qcov.txt')
  io.add_argument('--bed',dest = 'bed', help = 'PLINK2 binaries',
    default = '../params/bed')
  io.add_argument('--grm', dest = 'grm', help = 'Genetic correlation matrix',
    default = '/home/yh464/rds/rds-rb643-ukbiobank2/Data_Users/yh464/params/sp0.05_grm')
  
  params = parser.add_argument_group(title = 'parameters for GCTA')
  params.add_argument('--gcta', dest = 'gcta', help = 'Location of GCTA executable',
    default = '../toolbox/gcta/gcta64')
  params.add_argument('--maf', dest = 'maf', help = 'Filter by minor allele frequency',
    default = '0.01', type = str)
  params.add_argument('--keep', dest = 'keep', help = 'Subjects to keep', # intentionally absolute
    default = '/home/yh464/rds/rds-rb643-ukbiobank2/Data_Users/yh464/params/ukbkeepfile_202402.txt')
  params.add_argument('--extract', nargs='*', help = 'SNPs to extract from input files', default = [])
  
  xchr = parser.add_argument_group(title = 'X chromosome GWAS options')
  xchr.add_argument('--nox', dest = 'xchr', help = 'Do not conduct GWAS for X chromosome',
      default = True, action = 'store_false')
  xchr.add_argument('--xbed', help = 'PLINK binary for the X chromosome',
      default = '/home/yh464/rds/rds-rb643-ukbiobank2/Data_Users/yh464/params/bed/chrX')
  
  parser.add_argument('-f','--force', dest = 'force', help = 'Force overwrite',
    default = False, action = 'store_true')
  
  args = parser.parse_args()
  import os
  for arg in ['_in','out','gcta','dcov','qcov','grm','bed']:
      setattr(args, arg, os.path.realpath(getattr(args, arg)))
  
  from _utils import cmdhistory, path, logger
  logger.splash(args)
  cmdhistory.log()
  proj = path.project()
  proj.add_var('/%pheng',r'.+', 'phenotype group')
  proj.add_var('/%pheno',r'.+', 'phenotype')
  proj.add_var('/%maf',r'[0-9.]+', 'minor allele freq') # only allows digits and decimals
  proj.add_input(args._in+'/%pheng.txt', __file__)
  proj.add_output(args.out+'/%pheng/%pheno.fastGWA', __file__)
  try: main(args)
  except: cmdhistory.errlog()