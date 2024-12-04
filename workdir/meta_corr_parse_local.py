#!/usr/env/bin python
import argparse
parser = argparse.ArgumentParser(
  description = 'this script prepares the psychiatric rg for r-ggseg plotting')
parser.add_argument('pheno', nargs = '*', help = 'local phenotypes',
  default=['deg_local','degi_local','degc_local','clu_local','eff_local','mpl_local'])
parser.add_argument('-i','--in', dest = '_in', help = 'input directory',
  default='../local_corr/')
parser.add_argument('-o','--out', dest = 'out', help = 'output directory')
# always overwrites
args = parser.parse_args()

import os
args._in = os.path.realpath(args._in)
if type(args.out) == type(None): args.out = args._in

import pandas as pd

rg_list = []
z_list = []

os.chdir(args._in)
for x in args.pheno:
  df = pd.read_table(f'{x}/meta_rg_summary.csv',index_col = 0)
  df.insert(loc = 0, column = 'label', value = df.index)
  df = df.melt(id_vars = 'label', var_name = 'ext')
  df['pheno'] = x
  rg_list.append(df)
  
  df = pd.read_table(f'{x}/meta_z_summary.csv',index_col = 0)
  df.insert(loc = 0, column = 'label', value = df.index)
  df = df.melt(id_vars = 'label', var_name = 'ext')
  df['pheno'] = x
  z_list.append(df)

rg_df = pd.concat(rg_list,ignore_index = True)
for y in rg_df.ext.unique():
  tmp = rg_df.loc[rg_df['ext']==y, :].drop('ext', axis = 1)
  tmp = tmp.pivot(columns = 'pheno', index = 'label', values = 'value')
  tmp.columns.name = None
  tmp.to_csv(f'{args.out}/{y}_rg.csv')

z_df = pd.concat(z_list,ignore_index = True)
for y in z_df.ext.unique():
  tmp = z_df.loc[z_df['meta']==y, :].drop('meta', axis = 1)
  tmp = tmp.pivot(columns = 'pheno', index = 'label', values = 'value')
  tmp.columns.name = None
  tmp.to_csv(f'{args.out}/{y}_z.csv')