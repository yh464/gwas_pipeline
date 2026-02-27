#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def check_ancestry(reported, pc1, pc2, target):
  import pandas as pd
  filters = []
  for e in target:
    if e == 'eur':
      filters.append(reported.isin([1,1001,1002,1003]) & (pc1 > -39.574) & (pc1 < 16.714) & (pc2 > -17.145) & (pc2 < 23.645))
    elif e == 'mix':
      filters.append(reported.isin([2,2001,2002,2003,2004]) & (pc1 > -294.624) & (pc1 < 487.466) & (pc2 > -316.213) & (pc2 < 272.319))
    elif e == 'oas':
      filters.append(reported.isin([3,3001,3002,3003,3004]) & (pc1 > -36.806) & (pc1 < 191.373) & (pc2 > -304.746) & (pc2 < 74.524))
    elif e == 'afr':
      filters.append(reported.isin([4,4001,4002,4003]) & (pc1 > 79.454) & (pc1 < 649.381) & (pc2 > -47.610) & (pc2 < 180.314))
    elif e == 'eas':
      filters.append(reported.isin([5,5001,5002,5003]) & (pc1 > 95.421) & (pc1 < 215.474) & (pc2 > -377.089) & (pc2 < -159.840))
  filters = pd.concat(filters, axis = 1).any(axis = 1)
  return filters

def main(args):
  from time import perf_counter as t
  from fnmatch import fnmatch
  import os
  import pandas as pd
  from io import StringIO
  tic = t()
  if args.subj != 'all': 
      subj = open(args.subj,'r').read().splitlines()
      for i in range(len(subj)):
          if subj[i][:3] == 'UKB': subj[i] = subj[i][3:]
      subj = [int(s) for s in subj]
  else: subj = 'all'
  if len(args.pheno) > 0 and os.path.isfile(args.pheno[0]): pheno = open(args.pheno[0],'r').read().splitlines()
  else: pheno = args.pheno
  fout = open(args.out+'.txt','w')
  fin = open(args._in,'r')
  
  n = len(subj) if subj != 'all' else 'all'
  log.log(f'requesting data for {n} subjects and {len(pheno)} phenotypes')
  log.log(f'phenotypes to be extracted: {pheno}')
  
  # input data is in format something\tsomething\tsomething
  hdr = fin.readline().replace('\n','').split('\t')

  # initialise output columns, 
  valid_cols = ['FID','IID']
  valid_col_ids = [0,0] # eid is the first column in the UKB file header
  valid_col_match = [-1,-1] # this is for the log output
  
  # qc columns
  qc_cols = ['f.31.0.0', # sex, self-reported
      'f.21000.0.0', # ethnicity, self-reported
      'f.22001.0.0', # sex, genetic
      'f.22009.0.1','f.22009.0.2', # first two genetic PCs
      'f.22027.0.0'] # excessive heterozygosity
  qc_col_ids = []
  
  # filter header
  for i in range(len(hdr)):
    tmp = hdr[i]
    for j in range(len(pheno)):
      if fnmatch(tmp,'f.'+pheno[j]+'.*'):
        valid_cols.append(tmp)
        valid_col_ids.append(i)
        valid_col_match.append(j)
    for j in qc_cols:
      if fnmatch(tmp, j):
        qc_col_ids.append(i)
  
  # columns given in args.pheno not found in the UKB extract
  error_list = []
  for j in range(len(pheno)):
    if not (j in valid_col_match):
      error_list.append(j)
  if len(error_list) > 0:
    log.log('following phenotypes are not found in this UKB fetch:')
    for j in error_list:
      log.log(pheno[j])
    log.log('\n')
  
  print('\t'.join(valid_cols), file = fout)
  cols_to_read = sorted(set(valid_col_ids + qc_col_ids))
  cols_to_read_hdr = [hdr[i] for i in cols_to_read]
  toc = t() - tic
  
  # for quality control
  n_cols = len(hdr)
  log.log(f'processed headers. time = {toc:.2f} seconds')
  log.log(f'header contains {n_cols} columns')
  
  count = 0
  sex = 0
  eth = 0
  het = 0
  chunk = 0
  while True:
    # read 5000 lines at a time
    chunk += 1
    lines = StringIO()
    for i, line in enumerate(fin):
      lines.write(line)
      if (i+1) % 5000 == 0: break
    else: break # from the while loop only if the for loop is not broken, i.e. EOF is reached
    lines.seek(0)
    df = pd.read_table(lines, header = None, names = cols_to_read_hdr, usecols = cols_to_read, index_col = hdr[0])
    # first filter for subjects in query
    if subj != 'all':
      df = df.loc[df.index.isin(subj),:]
      subj = set(subj) - set(df.index)

    # then QC subjects
    if not args.noqc:
      sex_filter = df[qc_cols[0]] == df[qc_cols[2]] # genetic sex ~ self-reported
      het_filter = df[qc_cols[5]].isna() # heterozygosity, must be NA
      eth_filter = check_ancestry(df[qc_cols[1]], df[qc_cols[3]], df[qc_cols[4]], args.ethnicity)
      df = df.loc[sex_filter & het_filter & eth_filter,:]
      sex += (~sex_filter).sum()
      eth += (sex_filter & ~eth_filter).sum()
      het += (sex_filter & eth_filter & ~het_filter).sum()
    
    # write out to output file
    df = df[valid_cols[2:]] # only keep the phenotype columns
    df.insert(0, 'IID', df.index)
    df.index.name = 'FID'
    df.to_csv(fout, sep = '\t', header = False)
    count += df.shape[0]
    log.log(f'processed {chunk * 5000} subjects. time = {t() - tic:.2f} seconds')
  
  if subj != 'all' and len(subj) > 0:
    log.log('following subjects are not found in this UKB fetch:')
    for j in subj:
      log.log(j)
    log.log('\n')
  
  log.log(f'{sex} subjects excluded due to reported sex != genetic sex')
  log.log(f'{eth} subjects excluded based on ethnicity')
  log.log(f'{het} subjects excluded due to excessive heterozygosity')
  log.log(f'Resultant sample size is {count} subjects')
  fout.close()
  return

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(
    description = 'this script extracts selected subjects and columns from ukb')
  parser.add_argument('-s','--subj', dest = 'subj', help = 'subjects list',
    default = '/home/yh464/rds/rds-rb643-ukbiobank2/Data_Users/yh464/params/subjlist_full.txt')
  parser.add_argument('-p','--pheno', dest = 'pheno', nargs = '*', help = 'phenotypes list',
    default = [])
  parser.add_argument('-e', '--ethnicity', choices = ['eur','mix','afr','oas','eas'], default = ['eur'],nargs='*')
  parser.add_argument('--noqc', action = 'store_true', default = False, help = 'No QC at subject level')
  parser.add_argument('-i','--in', dest = '_in', help = 'input ukb fetch file, TAB format, NOT csv',
    default = '/home/yh464/rds/rds-rb643-ukbiobank2/Data_Phenotype/DataFetch_20022024/ukb677594.tab')
  parser.add_argument('-o','--out', dest = 'out', help = 'output prefix', required = True)
  args = parser.parse_args()
  import os
  for arg in ['_in','out']:
      setattr(args, arg, os.path.realpath(getattr(args, arg)))
  if args.subj != 'all': args.subj = os.path.realpath(args.subj)
  if len(args.pheno) > 0 and os.path.isfile(args.pheno[0]): args.pheno = os.path.realpath(args.pheno[0])

  from _utils import cmdhistory, logger
  cmdhistory.log()
  log = logger.logger(fname = f'{args.out}.log')
  try: main(args)
  except: cmdhistory.errlog()