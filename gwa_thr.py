#!/usr/bin/env python
import argparse
parser = argparse.ArgumentParser(
  description = 'This script filters out the most significant SNPs for regional traits')
parser.add_argument('-i','--in', dest = '_in', help = 'input GWA file')
parser.add_argument('-t','--thr', dest= 'thr',
  help = 'significance threshold. format: 5e-8,1e-5 (separate different thresholds w/ commas)',
  default = '3.1075e-11,5e-8')
parser.add_argument('-o','--out',dest = 'out', help = 'output directory')
parser.add_argument('-f','--force', dest = 'force', help = 'force overwrite',
                    action = 'store_true', default = False)
args = parser.parse_args()

import pandas as pd
import os

# parse input
prefix = args._in.split('/')[-1].replace('.fastGWA','')
thr = args.thr.split(',')
for i in range(len(thr)):
  thr[i] = float(thr[i])
thrmax = max(thr)

# stop if all output files exist
q = True
for x in thr:
  if not os.path.isfile(f'{args.out}/siglist.{x}.{prefix}.txt'): q = False

if q and (not args.force): quit()

# threshold GWA sumstats
df = pd.read_csv(args._in, sep = '\s+')
df_t = df.loc[df['P']<thrmax,:]
df_t.insert(loc = 0, column = 'node', value = prefix)

for x in thr:
  tmp = df_t.loc[df_t['P'] < x,:]
  out_fname = f'{args.out}/siglist.{x}.{prefix}.txt'
  if tmp.size == 0:
    open(out_fname,'w').close()
  else:
    tmp.to_csv(out_fname, sep = '\t', index = False, header = False)
  