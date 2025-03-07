#!/usr/bin/env python3
'''
Merges X and autosomes for single traits, and plots Manhattan plots
'''

import argparse
parser = argparse.ArgumentParser(description = 
  'This programme compiles Manhattan plots for all fastGWA output files for a single phenotype file')
parser.add_argument('pheno', help = 'Phenotype')
parser.add_argument('--file', dest = 'file', help = 'Input fastGWA file')
parser.add_argument('-i','--in', dest = '_in', help = 'GWA file directory',
  default = '../gwa/')
parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
  default = '../gwa/manhattan/')
parser.add_argument('-f','--force',dest = 'force', help = 'force output',
  default = False, action = 'store_true')
args = parser.parse_args()

import os
for arg in ['_in','out']:
    exec(f'args.{arg} = os.path.realpath(args.{arg})')

import time

tic = time.perf_counter()

import pandas as pd
from qmplot import manhattanplot, qqplot
import matplotlib.pyplot as plt

toc = time.perf_counter() - tic
print(f'Loaded modules. Time = {toc:.3f} seconds')

os.chdir(f'{args._in}/{args.pheno}')
x = args.file
out_fname = f'{args.out}/{args.pheno}/{x}'.replace('.fastGWA','.manhattan.pdf')
_,ax = plt.subplots(figsize = (6,2), constrained_layout = True)

if (not os.path.isfile(out_fname)) or args.force:
  df = pd.read_table(x)
  sig = df.loc[df.P < 5e-8,:]
  df_sig = []
  for chrom, pos in zip(sig.CHR, sig.POS):
      df_sig.append(df.loc[(df.CHR == chrom) & (df.POS > pos - 2e5) & 
                           (df.POS < pos + 2e5),:])
  df_0 = df.loc[df.P < 0.003,:]
  df_1 = df.loc[(df.P >= 0.001) & (df.P < 0.003),:]
  df_2 = df.loc[(df.P >= 0.003) & (df.P < 0.01),:]
  df_3 = df.loc[(df.P >= 0.01) & (df.P < 0.03),:]
  df_4 = df.loc[(df.P >= 0.03) & (df.P < 0.1),:]
  df_5 = df.loc[(df.P >= 0.1) & (df.P < 0.3),:]
  df_6 = df.loc[(df.P >= 0.3)]
  df = pd.concat([df_0,df_1.iloc[::int(df_1.shape[0]/4000),:],
                  df_2.iloc[::int(df_2.shape[0]/4000),:], 
                  df_3.iloc[::int(df_3.shape[0]/4000),:], 
                  df_4.iloc[::int(df_4.shape[0]/4000),:], 
                  df_5.iloc[::int(df_5.shape[0]/4000),:], 
                  df_6.iloc[::int(df_6.shape[0]/4000),:]] + df_sig).drop_duplicates()
  print(df.shape)
  df.sort_values(by = ['CHR','POS'], inplace = True)
  _,ax = plt.subplots(figsize = (6,2))
  manhattanplot(data = df,
                chrom = 'CHR',
                pos = 'POS',
                pv = 'P',
                snp = 'SNP',
                is_annotate_topsnp=True, # annotate sig. SNPs
                sign_marker_p = 5e-8,  # Genome wide significant p-value
                sign_marker_color="r",
                logp = True,
                ld_block_size = 1000000,
                text_kws = {'fontfamily': 'sans-serif', 'fontsize': 20},
                ax = ax)
  xtick = list(range(9)) + [10,12,14,17,20]
  ax.set_xticks(ax.get_xticks()[xtick], [x+1 for x in xtick])
  plt.savefig(out_fname, bbox_inches = 'tight')
  plt.savefig(out_fname.replace('.pdf','.png'), dpi = 400)
  plt.close()
  qqplot(data = df['P'], title = x.replace('.fastGWA',''),
         marker= '.', xlabel=r"Expected $-log_{10}{(P)}$",
           ylabel=r"Observed $-log_{10}{(P)}$")
  plt.savefig(out_fname.replace('.manhattan','.qqplot'), bbox_inches = 'tight')
  plt.close()
  toc = time.perf_counter()-tic
  print(f'Fig plotted, time = {toc:.3f} seconds.')