#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2024-11-12

Conducts MAGMA for single GWAS summary statistics
To be used in combination with annot_magma_batch.py

Requires following inputs: 
    GWAS summary statistics (single file)
    MAGMA annotation files (including H-MAGMA)
    MAGMA binary
'''

import argparse
parser = argparse.ArgumentParser(
  description = 'This file constitutes the MAGMA pipeline for a single IDP')
parser.add_argument('pheno', help = 'Phenotype')
parser.add_argument('-i','--in', dest = '_in', help = 'Input GWA summary statistic file')
parser.add_argument('-d','--dir', dest = 'dir', help = 'Directory containing all phenotypes',
  default = '../gwa/')
parser.add_argument('-a', '--annot', dest ='annot', help = 'directory to annotation files',
  default = '../toolbox/hmagma')
parser.add_argument('-b', '--bfile', dest = 'bfile', help = 'bed binary to use in magma',
  default = '../toolbox/magma/g1000_eur')
parser.add_argument('-o', '--out', dest = 'out', help = 'output directory',
  default = '../annot/magma/')
parser.add_argument('--magma', dest = 'magma', help = 'MAGMA executable',
  default = '../toolbox/magma/magma')
parser.add_argument('--gset', dest ='gset', help = 'Gene sets to study enrichment',
  default = '../params/gset.txt')
parser.add_argument('-f','--force',dest = 'force', help = 'force output',
  default = False, action = 'store_true')
args = parser.parse_args()

import os
from fnmatch import fnmatch
# path normalisation
for arg in ['dir','annot','bfile','out','magma','gset']:
    exec(f'args.{arg} = os.path.realpath(args.{arg})')

os.system('module load gcc/11') # IMPORTANT!!!!!
os.chdir(f'{args.dir}/{args.pheno}')
prefix = args._in.replace('.fastGWA','')

annot_list = []
for x in os.listdir(args.annot):
    if fnmatch(x, '*.genes.annot'): 
        annot_list.append(x)

genes_dir = f'{args.out}/{args.pheno}/genes'
gsa_dir = f'{args.out}/{args.pheno}/gsa'
if not os.path.isdir(genes_dir): os.system(f'mkdir -p {genes_dir}')
if not os.path.isdir(gsa_dir): os.system(f'mkdir -p {gsa_dir}')

for annot in annot_list:
    # magma gene annotation
    out_geneset = f'{genes_dir}/{prefix}_{annot}'
    if (not os.path.isfile(f'{out_geneset}.genes.raw')) or args.force:
      os.system(f'{args.magma} --bfile {args.bfile} --gene-annot {args.annot}/{annot} '+
        f'--pval {args._in} ncol=N --out {out_geneset}')
      
    f = open(args.gset)
    flist = f.read().splitlines()
    
    # magma enrichment analysis
    for x in flist:
      xprefix = x.split('/')[-1].replace('.txt','')
      gsaout = f'{gsa_dir}/{prefix}.{annot}.{xprefix}'
      
      if (not os.path.isfile(f'{gsaout}.gsa.out')) or args.force:
        os.system(f'{args.magma} --gene-results {out_geneset}.genes.raw --set-annot'+
                  f' {x} --out {gsaout}')