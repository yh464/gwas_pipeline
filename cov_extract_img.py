#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-05-27
Extracts imaging covariates given a list of subjects
'''
def parse_euler(eid):
  import os; import numpy as np
  euler_file = f'{eid}/surfaces/{eid}/stats/aseg.stats'
  if not os.path.isfile(euler_file): return np.nan
  lines = open(euler_file).read().splitlines()
  for line in lines:
    if line.startswith('# Measure SurfaceHoles'): return int(line.split(', ')[-2])
  return np.nan

def parse_fd(eid):
  import os; import numpy as np
  fd_file = f'{eid}/func/fMRI/rfMRI.ica/mc/prefiltered_func_data_mcf.par'
  if not os.path.isfile(fd_file): return np.nan, np.nan
  df = np.loadtxt(fd_file)
  df[:,:3] *= np.pi * 50 / 180 # rotations
  fd = df[2:,:] - df[:-2,:]
  fd = np.abs(fd).sum(axis = 1)
  return np.sum(fd)/(len(fd) + 1), np.max(fd)

def main(args):
  import pandas as pd
  subj = open(args.subj).read().splitlines()
  subj = [x.split()[0] for x in subj]
  os.chdir(args._in)
  covs = []
  for x in subj:
    if not x.startswith('UKB'): x = f'UKB{x}'
    fd, fdmax = parse_fd(x)
    covs.append(pd.DataFrame(dict(
      FID = [x.replace('UKB','')], IID = x.replace('UKB',''),
      euler = parse_euler(x),
      fd = fd, fdmax = fdmax
    )))
  covs = pd.concat(covs, ignore_index = True).dropna()
  covs.to_csv(f'{args.out}.txt', sep = '\t', index = False, header = True)

if __name__ == '__main__':
  from _utils.slurm import slurm_parser
  parser = slurm_parser(description = 'Extracts imaging covariates given a list of subjects')
  parser.add_argument('-i','--in', dest = '_in', default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Imaging',
    help = 'Root directory of imaging files, contains one subfolder per subject')
  parser.add_argument('-s','--subj', help = 'subjects list, one per line or BIM format',
    default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/ukbkeepfile_202002.txt')
  parser.add_argument('-o','--out', help = 'Output prefix')
  # always overwrites
  args = parser.parse_args()
  import os
  args._in = os.path.abspath(args._in); args.out = os.path.abspath(args.out)
  from _utils import cmdhistory, logger
  logger.splash(args)
  cmdhistory.log()
  try: main(args)
  except: cmdhistory.errlog()
  