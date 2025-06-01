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

def main(args):
  # array submitter
  from _utils.slurm import array_submitter
  submitter = array_submitter(name = f'gwa_{args.pheno}',timeout = 90)
  
  # locate phenotype file
  import os
  import fnmatch
  flist = []
  for f in os.listdir(args._in):
    if fnmatch.fnmatch(f,f'*{args.pheno}*') and not(os.path.isdir(f)):       # search for all files matching args.pheno
      flist.append(f'{args._in}/{f}')
  if len(flist) != 1: raise ValueError('Please give only ONE phenotype file')
  
  # temp and log
  tmpdir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/temp/'                                 # temporatory dir
  if not os.path.isdir(tmpdir): os.mkdir(tmpdir)
  
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

  # check validity of the phenotype file
  import pandas as pd
  f = flist[0]
  if not fnmatch.fnmatch(f,'*.txt'): raise ValueError('Phenotype file should be in TXT format')
  os.system(f'head {f} -n 5 > {tmpdir}/temp_{args.pheno}.txt')
  df = pd.read_csv(f'{tmpdir}/temp_{args.pheno}.txt',sep = '\s+')
  c = df.columns.values
  if c[0] != 'FID' or c[1] != 'IID':
    raise ValueError('Phenotype file should be in the format: FID IID *pheno')
  os.remove(f'{tmpdir}/temp_{args.pheno}.txt')
  
  # create output folder
  outdir = f'{args.out}/{os.path.basename(f)}/'.replace('.txt','')
  print(outdir)
  if not os.path.isdir(outdir):
    os.system(f'mkdir -p {outdir}')                                              # this also generates args.out
  
  # phenotypes to be analysed
  c = c[2:]
  print('Following traits are to be GWA-analysed:')
  for i in c: print(i)
  
  # for each phenotype
  for i in range(c.size):
    mpheno = i+1
    trait = c[i]
    out_fname = outdir + trait
    # check existing files
    if os.path.isfile(f'{out_fname}.fastGWA') and not args.force:
      print(f'Trait already analysed for: {trait}')
      continue
    
    submitter.add(
      f'python gwa_by_trait.py -i {f} -o {out_fname} --mpheno {mpheno} --dcov {args.dcov} '+
      f'--qcov {args.qcov} --bed {args.bed} --grm {args.grm} --gcta {args.gcta} --maf {args.maf} '+
      f'--keep {args.keep} {xchr} --xbed {args.xbed} {extract} {force}'
      )
  submitter.submit()

if __name__ == '__main__':
  import argparse
  from _utils.slurm import parser_config
  # argument input
  parser = argparse.ArgumentParser(description=
    'This programme runs GWA for any phenotype given as the 1st positional argument')
  parser.add_argument('pheno', help = 'Phenotype file in TXT format - please supply ONLY ONE')
  
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
    default = '../params/bed_files_ukb.txt')
  io.add_argument('--grm', dest = 'grm', help = 'Genetic correlation matrix',
    default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/sp0.05_grm')
  
  params = parser.add_argument_group(title = 'parameters for GCTA')
  params.add_argument('--gcta', dest = 'gcta', help = 'Location of GCTA executable',
    default = '../toolbox/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1')
  params.add_argument('--maf', dest = 'maf', help = 'Filter by minor allele frequency',
    default = '0.01', type = str)
  params.add_argument('--keep', dest = 'keep', help = 'Subjects to keep', # intentionally absolute
    default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/ukbkeepfile_202402.txt')
  params.add_argument('--extract', nargs='*', help = 'SNPs to extract from input files', default = [])
  
  xchr = parser.add_argument_group(title = 'X chromosome GWAS options')
  xchr.add_argument('--nox', dest = 'xchr', help = 'Do not conduct GWAS for X chromosome',
      default = True, action = 'store_false')
  xchr.add_argument('--xbed', help = 'PLINK binary for the X chromosome',
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/bed/chrX')
  
  parser.add_argument('-f','--force', dest = 'force', help = 'Force overwrite',
    default = False, action = 'store_true')
  
  args = parser.parse_args()
  import os
  for arg in ['_in','out','gcta','dcov','qcov','grm','bed']:
      exec(f'args.{arg} = os.path.realpath(args.{arg})')
  
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