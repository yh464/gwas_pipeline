#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def main(args):
  import os
  from fnmatch import fnmatch
  import pandas as pd
  import seaborn as sns
  import matplotlib.pyplot as plt
  
  os.chdir(args._in)
  if not os.path.isdir('diagnostics'): os.mkdir('diagnostics')
  for f in os.listdir():
    if not fnmatch(f, '*.txt'): continue
    prefix = os.path.basename(f).replace('.txt','')
    out_fname = f'diagnostics/{prefix}.png'
    if os.path.isfile(out_fname) and not args.force: continue
    df = pd.read_csv(f, sep = args.sep).melt(id_vars = ['FID','IID'])
    sns.displot(data = df, x = 'value', col = 'variable', height = 10)
    plt.savefig(out_fname)
    plt.close()
    print(f'Processed: {f}')

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description = 'This script batch plots histograms for phenotypes')
  parser.add_argument('-i','--in', dest = '_in', help = 'Input directory', 
    default = '../pheno/ukb/')
  parser.add_argument('-s','--sep', dest = 'sep', help = 'column separator', 
    default = '\s+')
  parser.add_argument('-f','--force', dest = 'force', help = 'force overwrite',
    action = 'store_true', default = False)
  args = parser.parse_args()
  import os
  args._in = os.path.realpath(args._in)
  
  from ._utils import cmdhistory, path
  cmdhistory.log()
  proj = path.project()
  proj.add_input(args._in, __file__)
  proj.add_output(args._in+'/diagnostics/', __file__)
  try: main(args)
  except: cmdhistory.errlog()