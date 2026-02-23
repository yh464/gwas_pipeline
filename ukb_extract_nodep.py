#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def check_ancestry(reported, pc1, pc2, target):
  for e in target:
    if e == 'eur' and reported in ['1','1001','1002','1003'] and \
      -39.574 < pc1 < 16.714 and -17.145 < pc2 < 23.645: return True
    elif e == 'mix' and reported in ['2','2001','2002','2003','2004'] and \
      -294.624 < pc1 < 487.466 and -316.213 < pc2 < 272.319: return True
    elif e == 'oas' and reported in ['3','3001','3002','3003','3004'] and \
      -36.806 < pc1 < 191.373 and -304.746 < pc2 < 74.524: return True
    elif e == 'afr' and reported in ['4','4001','4002','4003'] and \
      79.454 < pc1 < 649.381 and -47.610 < pc2 < 180.314: return True
    elif e == 'eas' and reported == '5' and \
      95.421 < pc1 < 215.474 and -377.089 < pc2 < -159.840: return True
  return False

def main(args):
  from time import perf_counter as t
  from fnmatch import fnmatch
  import os
  tic = t()
  idx = 0
  if args.subj != 'all': 
      subj = open(args.subj,'r').read().splitlines()
      for i in range(len(subj)):
          if subj[i][:3] == 'UKB': subj[i] = subj[i][3:]
  else: subj = 'all'
  if len(args.pheno) > 0 and os.path.isfile(args.pheno[0]): pheno = open(args.pheno[0],'r').read().splitlines()
  else: pheno = args.pheno
  fout = open(args.out+'.txt','w')
  flog = open(args.out+'.log','w')
  fin = open(args._in,'r')
  
  n = len(subj) if subj != 'all' else 'all'
  print(f'requesting data for {n} subjects and {len(pheno)} phenotypes\n')
  print(f'phenotypes to be extracted: {pheno}')
  print(f'requesting data for {n} subjects and {len(pheno)} phenotypes\n', file = flog)
  print(f'phenotypes to be extracted: {pheno}', file = flog)
  
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
    print('following phenotypes are not found in this UKB fetch:', file = flog)
    for j in error_list:
      print(pheno[j], file = flog)
    print('\n', file = flog)
  
  print('\t'.join(valid_cols), file = fout)
  toc = t() - tic
  
  # for quality control
  n_cols = len(hdr)
  print(f'processed headers. time = {toc:.2f} seconds')
  print(f'processed headers. time = {toc:.2f} seconds', file = flog)
  print(f'header contains {n_cols} columns')
  print(f'header contains {n_cols} columns', file = flog)
  
  count = 0
  sex = 0
  eth = 0
  het = 0
  while True:
    line = fin.readline()
    if len(line) == 0: break
    idx += 1
    if idx % 5000 == 0:
      toc = t() - tic
      print(f'read {idx} subjects. time = {toc:.2f} seconds')
      print(f'read {idx} subjects. time = {toc:.2f} seconds', file = flog)
      
    # split columns
    line = line.replace('\n','').split('\t')
    
    if len(line) != n_cols:
        print(f'{line[0]}: length does not match header, {len(line)} columns')
        print(f'{line[0]}: length does not match header, {len(line)} columns', file = flog)
        continue
    if subj != 'all' and not line[0] in subj: continue # skip unwanted subjects
    if subj != 'all': subj.remove(line[0])
    
    if not args.noqc:
      # first QC subjects
      if not line[qc_col_ids[0]] == line[qc_col_ids[2]]: sex += 1; continue # genetic sex ~ self-reported
      if not fnmatch(line[qc_col_ids[5]],'[Nn][Aa]'): het += 1; continue # heterozygosity, must be NA
      
      if line[qc_col_ids[3]] == 'NA': continue
      pc1 = float(line[qc_col_ids[3]])
      if line[qc_col_ids[4]] == 'NA': continue
      pc2 = float(line[qc_col_ids[4]])
      if not check_ancestry(line[qc_col_ids[1]], pc1, pc2, args.ethnicity): eth += 1; continue
    # if QC is passed, then write out to output file
    line_out = [line[i] for i in valid_col_ids]
    print('\t'.join(line_out), file = fout)
    count += 1
  
  if subj != 'all' and len(subj) > 0:
    print('following subjects are not found in this UKB fetch:', file = flog)
    for j in subj:
      print(j, file = flog)
    print('\n', file = flog)
  
  print(f'{sex} subjects excluded due to reported sex != genetic sex')
  print(f'{eth} subjects excluded based on ethnicity')
  print(f'{het} subjects excluded due to excessive heterozygosity')
  print(f'Resultant sample size is {count} subjects')
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
      exec(f'args.{arg} = os.path.realpath(args.{arg})')
  if args.subj != 'all': args.subj = os.path.realpath(args.subj)
  if len(args.pheno) > 0 and os.path.isfile(args.pheno[0]): args.pheno = os.path.realpath(args.pheno[0])
  
  main(args)
