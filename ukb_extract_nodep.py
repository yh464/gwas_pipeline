#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def main(args):
  from time import perf_counter as t
  from fnmatch import fnmatch
  import os
  tic = t()
  idx = 0
  subj = open(args.subj,'r').read().splitlines()
  if args.pheno == 'null': pheno = [] # just quality control
  elif os.path.isfile(args.pheno): pheno = open(args.pheno,'r').read().splitlines()
  else: pheno = args.pheno.split(',')
  fout = open(args.out+'.txt','w')
  flog = open(args.out+'.log','w')
  fin = open(args._in,'r')
  
  print(f'requesting data for {len(subj)} subjects and {len(pheno)} phenotypes=\n', file = flog)
  print(f'phenotypes to be extracted: {pheno}', file = flog)
  
  # process subjects list
  for i in range(len(subj)):
      if subj[i][:3] == 'UKB': subj[i] = subj[i][3:]
  
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
      if fnmatch(tmp,'f.'+pheno[j]+'*'):
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
    # line[0] = line[0].replace('"','')
    # line[-1] = line[-1].replace('"\n','')
    
    # some lines contain a sub-table in csv format, so this segment screens for those
    # sub-tables contain the character '"' which should not appear elsewhere
    # for tmp in range(len(line)-1, 0, -1): # reverse order because list length will change
    #     if line[tmp-1][:2] == ',"':
    #         line[tmp-1] = line[tmp-1][2:]
    #         line[tmp-2] += ',"' # this is the end of a sub-table
    
    #     if line[tmp][-2:] == ',"' and line[tmp-1].find('"') == 0 \
    #         and line[tmp-1][-2:] != ',"':
    #         line[tmp] = '","'.join(line[tmp-1:tmp+1])
    #         del line[tmp-1]
    
    if len(line) != n_cols:
        print(f'{line[0]}: length does not match header, {len(line)} columns')
        print(f'{line[0]}: length does not match header, {len(line)} columns', file = flog)
        continue
    if not line[0] in subj: continue # skip unwanted subjects
    subj.remove(line[0])
    
    # first QC subjects
    if not line[qc_col_ids[1]] in ['1','1001','1002','1003','1004']: continue # European ancestry
    if not line[qc_col_ids[0]] == line[qc_col_ids[2]]: continue # genetic sex ~ self-reported
    if not fnmatch(line[qc_col_ids[5]],'[Nn][Aa]'): continue # heterozygosity, must be NA
    
    if line[qc_col_ids[3]] == 'NA': continue
    pc1 = float(line[qc_col_ids[3]])
    if pc1 < -272.541 or pc1 > 270.278: continue # first genetic PC, 5SD is manually calculated
    
    if line[qc_col_ids[4]] == 'NA': continue
    pc2 = float(line[qc_col_ids[4]])
    if pc2 < -138.931 or pc2 > 139.650: continue # second genetic PC, 5SD manually calculated
    
    # if QC is passed, then write out to output file
    line_out = [line[i] for i in valid_col_ids]
    print('\t'.join(line_out), file = fout)
  
  if len(subj) > 0:
    print('following subjects are not found in this UKB fetch:', file = flog)
    for j in subj:
      print(j, file = flog)
    print('\n', file = flog)
  
  return

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(
    description = 'this script extracts selected subjects and columns from ukb')
  parser.add_argument('-s','--subj', dest = 'subj', help = 'subjects list',
    default = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/params/subjlist_rsfmri_hcp.txt')
  parser.add_argument('-p','--pheno', dest = 'pheno', help = 'phenotypes list, file or string split by ","',
    default = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/params/ukb_pheno.txt')
  parser.add_argument('-i','--in', dest = '_in', help = 'input ukb fetch file, TAB format, NOT csv',
    default = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Phenotype/DataFetch_20022024/ukb677594.tab')
  parser.add_argument('-o','--out', dest = 'out', help = 'output prefix', required = True)
  args = parser.parse_args()
  
  main(args)
