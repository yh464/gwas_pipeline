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
parser.add_argument('-p','--plot-only', dest = 'p', help = 'skip concatenation, only plot graphs',
  default = False, action = 'store_true')
parser.add_argument('-a','--autosome-only',dest = 'a', help = 'exclude sex chromosomes',
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
out_fname = f'{args.out}/{args.pheno}/{x}'.replace('.fastGWA','.manhattan.png')
_,ax = plt.subplots(figsize = (36,12), constrained_layout = True)

if ((not os.path.isfile(x.replace('.fastGWA','_all_chrs.fastGWA'))) or args.force) \
  and (not args.a):
  df_a = pd.read_csv(x,sep = '\t')
  df_x = pd.read_csv(x.replace('.fastGWA','_X.fastGWA'), sep = '\t')
  df = pd.concat([df_a, df_x], axis = 0).sort_values(by = ['CHR','POS'])
  df.to_csv(x.replace('.fastGWA','_all_chrs.fastGWA'), sep = '\t', index = False)
toc = time.perf_counter() - tic
print(f'GWA files concatenated, time = {toc:.3f} seconds')

if (not os.path.isfile(out_fname)) or args.force or args.p:
  if args.a:
    df = pd.read_csv(x, sep = '\t')
  else:
    df = pd.read_csv(x.replace('.fastGWA','_all_chrs.fastGWA'), sep = '\t')
  df.sort_values(by = ['CHR','POS'], inplace = True)
  _,ax = plt.subplots(figsize = (12,4))
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
                text_kws = {'fontfamily': 'sans-serif', 'fontsize': 10},
                ax = ax)
  # fig = plt.gcf()
  # fig.set_size_inches(12,4)
  # plt.title(x.replace('.fastGWA','').replace('_0.01',''),fontsize = 14)
  plt.savefig(out_fname, dpi = 400)
  plt.close()
  qqplot(data = df['P'], title = x.replace('.fastGWA',''),
         marker= '.', xlabel=r"Expected $-log_{10}{(P)}$",
           ylabel=r"Observed $-log_{10}{(P)}$")
  plt.savefig(out_fname.replace('manhattan.png','qqplot.png'), dpi = 400)
  plt.close()
  toc = time.perf_counter()-tic
  print(f'Fig plotted, time = {toc:.3f} seconds.')