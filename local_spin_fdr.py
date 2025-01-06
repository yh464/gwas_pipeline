#!/usr/bin/env python
# -*- coding: utf-8 -*-

def main(args):
  import os
  import pandas as pd
  from scipy.stats import false_discovery_control as fdr
  from fnmatch import fnmatch
  
  os.chdir(args._in)
  for x in os.listdir():
    if not fnmatch(x, '*.csv'): continue
    df = pd.read_csv(x, sep = '\s+')
    try:
      p = df['p1tail']
    except: continue
    q = fdr(p)
    df['fdr1tail'] = q
    df['Class'] = df['Class'].replace('7Networks_1', 'Visual').replace(
      '7Networks_2', 'Somatomotor').replace('7Networks_3','D_Attention').replace(
      '7Networks_4','V_Attention').replace('7Networks_5','Limbic').replace(
      '7Networks_6','Frontoparietal').replace('7Networks_7','Default')
    df.to_csv(x, sep = '\t', index = False)

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description = 'FDR script for local spin perm tests')
  parser.add_argument('-i','--in', dest = '_in', help = 'input directory',
    default = '../local_corr/')
  args = parser.parse_args()
  import os
  args._in = os.path.realpath(args._in)
  
  from _utils import cmdhistory, path
  cmdhistory.log()
  proj = path.project()
  proj.add_input(args._in+'/%pheno_.*.csv', __file__)
  proj.add_output(args._in+'/%pheno_.*.csv', __file__) # overwrites file!
  try: main(args)
  except: cmdhistory.errlog()