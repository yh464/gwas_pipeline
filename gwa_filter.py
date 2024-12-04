#!/usr/bin/env python3
'''
This scripts filters fastGWA files
'''

import argparse

# argument input
parser = argparse.ArgumentParser(description=
  'This programme filters the GRM to different thresholds')
parser.add_argument('-i', dest = '_in', help = 'input fastGWA file')
parser.add_argument('-o', dest = 'out', help = 'output fastGWA file')
parser.add_argument('--freq', dest = 'freq', help = 'MAF filter')
args = parser.parse_args()

import pandas as pd
df = pd.read_csv(args._in,sep = '\t')
df = df[df.AF1 >= float(args.freq)]
df = df[df.AF1 <= 1-float(args.freq)]
df.to_csv(args.out, index = False, sep = '\t')