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

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams['font.sans-serif'] = 'Arial'

toc = time.perf_counter() - tic
print(f'Loaded modules. Time = {toc:.3f} seconds')

os.chdir(f'{args._in}/{args.pheno}')
x = args.file
if not os.path.isdir(f'{args.out}/{args.pheno}'): os.mkdir(f'{args.out}/{args.pheno}')
out_fname = f'{args.out}/{args.pheno}/{x}'.replace('.fastGWA','.manhattan.pdf')
_,ax = plt.subplots(figsize = (6,2), constrained_layout = True)

if (not os.path.isfile(out_fname)) or args.force:
  df = pd.read_table(x).sort_values(by = ['CHR','POS'])
  sig = df.loc[df.P < 1e-3,:]
  
  # truncated Manhattan plot, pdf
  _,ax = plt.subplots(figsize = (6,2))
  manhattanplot(data = sig,
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
  plt.close()
  toc = time.perf_counter()-tic
  
  # full Manhattan plot
  _,ax = plt.subplots(figsize = (6,2))
  manhattanplot(data = sig,
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
  plt.savefig(out_fname.replace('pdf','png'), dpi = 500, bbox_inches = 'tight')
  plt.close()
  
  plt.rcParams['font.size'] = 20
  _, ax = plt.subplots(figsize = (3,3))
  qqplot(data = df['P'], title = '', ax = ax,
         marker= '.', xlabel=r"Expected $-log_{10}{(P)}$",
           ylabel=r"Observed $-log_{10}{(P)}$")
  plt.savefig(out_fname.replace('.manhattan.pdf','.qqplot.png'), dpi = 500, bbox_inches = 'tight')
  plt.close()
  print(f'Fig plotted, time = {toc:.3f} seconds.')