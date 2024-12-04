#!/usr/bin/env python3
'''
This script runs gene-based annotation for each finemap
'''

import argparse
parser = argparse.ArgumentParser(
  description = 'This script runs gene-based annotation for each finemap')
parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
parser.add_argument('-i','--in', dest = '_in', help = 'Directory containing all phenotypes', # in and out to same directory
  default = '../finemap/')
parser.add_argument('-o','--out', dest = 'out', help = 'Output directory') # in and out to same directory
parser.add_argument('--magma', dest = 'magma', help = 'MAGMA executable',
  default = '../toolbox/magma/magma')
parser.add_argument('-b', '--bfile', dest = 'bfile', help = 'bed binary to use in magma',
  default = '../toolbox/magma/g1000_eur')
parser.add_argument('-a', '--annot', dest ='annot', help = 'directory of annotation files',
  default = '../toolbox/hmagma')
parser.add_argument('-r', '--ref', dest = 'ref', help = 'Gene label document',
  default = '../toolbox/magma/ENSG.gene.loc')
parser.add_argument('-f','--force',dest = 'force', help = 'force overwrite',
  default = False, action = 'store_true')
args = parser.parse_args()
import os
for arg in ['_in','magma','bfile','annot','ref']:
    exec(f'args.{arg} = os.path.realpath(args.{arg})')
if type(args.out) == type(None): args.out = args._in
else: args.out = os.path.realpath(args.out)

# make output directories
from fnmatch import fnmatch
import pandas as pd

if not os.path.isdir(args.out): os.mkdir(args.out)
logdir = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/logs/'
log = open(f'{logdir}finemap_annot.log','w')
tmpdir = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/temp/finemap_annot_cache/'
if not os.path.isdir(tmpdir): os.mkdir(tmpdir)

# gene labelling
ref = pd.read_csv(args.ref, sep = '\t', header = None)
ref.columns = ['GENE','CHR','START','STOP','STRAND','LABEL']
ref_lbl = ref[['GENE','LABEL']]

annot_list = []
for x in os.listdir(args.annot):
    if fnmatch(x, '*.genes.annot'): annot_list.append(x)

for x in args.pheno:
  os.chdir(args._in)
  os.chdir(x)
  if not os.path.isdir('annot'): os.mkdir('annot')
  
  # prefix
  pflist = []
  for y in os.listdir():
    if fnmatch(y, '*.finemap.summary'):
      pflist.append(y.replace('.finemap.summary',''))
  
  os.chdir('polyfun_stats')
  all_phenos = []
  for y in pflist:
    print(f'Processing phenotype: {y}', file = log)
    print(f'Processing phenotype: {y}')
    all_gsets = []
    for annot in annot_list:
        # output filter
        out_fname = f'{args.out}/{x}/annot/{y}.{annot}.finemap.annot'
        if args.force:
          try: os.remove(out_fname)
          except: pass
        if os.path.isfile(out_fname):
          df = pd.read_csv(out_fname, sep = '\s+')
          df = df.merge(ref_lbl, how = 'inner', on = 'GENE')
          df['METHOD'] = annot
          all_gsets.append(df)
          df.to_csv(f'{out_fname}.lab', index = False)
          continue
        
        # all finemap outputs
        flist = []
        for z in os.listdir():
          if fnmatch(z, f'{y}*.csv'):
            flist.append(z)
        
        sig_sets = []
        for z in flist:
          print(f'>>>>> Processing segment: {z}', file = log)
          # 95% credible set
          zprefix = z.replace('.csv','')
          df = pd.read_csv(z, sep = '\s+')
          
          n_set = df['CREDIBLE_SET'].values.max()
          print(f'>>>>> {n_set} significant sets found for segment {z}', file = log)
          sig_sets.append(n_set)
          
          for i in range(1,n_set+1):
            df_out = df.loc[df.CREDIBLE_SET == i,:]
            df_out.to_csv(f'{tmpdir}/{zprefix}.sigset{i}.sigsnp', index = False, sep = '\t')
            
            # annotate by MAGMA
            out_geneset = f'{tmpdir}/{zprefix}.sigset{i}'
            o = 0
            
            if (not os.path.isfile(f'{out_geneset}.genes.raw')) or args.force:
              o = os.system(f'{args.magma} --bfile {args.bfile} --gene-annot {args.annot}/{annot} '+
                f'--pval {tmpdir}/{zprefix}.sigset{i}.sigsnp ncol=N --out {out_geneset}')
            
            if o != 0: 
              print(f'WARNING: magma failed for {out_geneset}')
              print(f'WARNING: magma failed for {out_geneset}', file = log)
              continue
            
            # concatenate output
            if not os.path.isfile(out_fname):
              os.system(f'cat {out_geneset}.genes.out > {out_fname}')
            else:
              os.system(f'tail -n +2 {out_geneset}.genes.out >> {out_fname}')
        
        if not os.path.isfile(out_fname): continue # if no result
        df = pd.read_csv(out_fname, sep = '\s+')
        df = df.merge(ref_lbl, how = 'inner', on = 'GENE')
        df.to_csv(f'{out_fname}.lab', index = False, sep = '\t')
        df['METHOD'] = annot.replace('.genes.annot','')
        all_gsets.append(df)
    
    # merge different annotations
    if len(all_gsets) == 0:
        print(f'WARNING: phenotype {y} does not have significant SNPs')
        continue
    
    all_gsets = pd.concat(all_gsets).drop_duplicates()
    all_gsets.to_csv(f'{args.out}/{x}/{y}.finemap.annot', sep ='\t', index = False)
    all_gsets.insert(0, 'PHENO',y)
    all_phenos.append(all_gsets)
    del all_gsets
  
  all_phenos = pd.concat(all_phenos).drop_duplicates().sort_values(by = ['LABEL','P'])
  all_phenos.to_csv(f'{args.out}/{x}_finemap_genes.txt', sep = '\t', index = False)
# clear temp
os.chdir(tmpdir)
# for x in os.listdir(): os.remove(x)