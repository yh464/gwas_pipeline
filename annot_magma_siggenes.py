#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2023-07-20
Version 2: 2024-11-14

Summarises gene level summary statistics for HMAGMA and MAGMA outputs

Preceding workflow:
    annot_batch.py
Requires following inputs:
    MAGMA genes.out outputs
Changelog:
    Added the 'annot' paramter to reflect the multitude of Hi-C and nearest-gene-based annotations
    Changed output to long-format table
'''

def main(args):
    import os
    from fnmatch import fnmatch
    import pandas as pd
    import scipy.stats as sts
    from _utils.path import normaliser
    norm = normaliser()
    
    # reference gene labelling
    ref = pd.read_table(args.ref)
    ref_lbl = ref[['GENE','LABEL']]
    
    # annotation method
    annot_list = []
    for x in os.listdir(args.annot):
        if fnmatch(x, '*.genes.annot'): 
            annot_list.append(x.replace('.genes.annot',''))
    
    for x in args.pheno:
      print(f'Processing: {x}')
      os.chdir(f'{args._in}/{x}')
      all_annot = []
      
      if not os.path.isdir(f'{args._in}/{x}/summary'):
          os.mkdir(f'{args._in}/{x}/summary')
      
      for annot in annot_list:
          dflist = []
          for y in os.listdir():
            if fnmatch(y, f'*{annot}.genes.out'):
              prefix = y.replace(f'_{annot}.genes.out','')
              df = pd.read_table(y, sep = '\\s+')
              df.drop(['CHR','START','STOP','NSNPS','NPARAM','N'], axis = 1, inplace = True)
              df = df.merge(ref_lbl, how = 'inner', on = 'GENE')
              tmp = df.P.values.copy()    
              fdr = sts.false_discovery_control(tmp)
              df['Pfdr']= fdr
              df.insert(0, 'annot', value = annot)
              df.insert(0, 'pheno', value = prefix)        
              dflist.append(df)
          
          summary = pd.concat(dflist).sort_values(by = 'Pfdr')
          summary.to_csv(f'{args._in}/{x}/summary/{annot}.siggenes.txt', sep = '\t', index = False)
          all_annot.append(summary)
      all_annot = pd.concat(all_annot, ignore_index = True)
      all_annot = all_annot.loc[all_annot.P <= 0.05 ,:].sort_values(by = ['Pfdr','P'])
      norm.normalise(all_annot).to_csv(f'{args._in}/{x}/all_siggenes.txt', sep = '\t', index = False)
      
      all_annot = all_annot.loc[all_annot.Pfdr <= 0.05,:]
      overlaps = pd.DataFrame(data = 0, index = all_annot.annot.unique(), columns = all_annot.LABEL.unique())
      for i in all_annot.index:
          overlaps.loc[all_annot.loc[i, 'annot'], all_annot.loc[i,'LABEL']] = 1
      overlaps.index.name = 'annot'
      norm.normalise(overlaps).to_csv(f'{args._in}/{x}/all_siggenes_overlaps.txt', sep = '\t', index = True, header = True)
      
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'Parses the raw ***.genes.out MAGMA output for each pheno')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
    parser.add_argument('-i','--in', dest = '_in', help = 'Input directory',
      default = '../annot/magma/')
    parser.add_argument('-a', '--annot', dest ='annot', help = 'directory to annotation files',
      default = '../toolbox/hmagma')
    parser.add_argument('-r', '--ref', dest = 'ref', help = 'Gene label document',
      default = '../params/genes_ref.txt')
    parser.add_argument('-f','--force',dest = 'force', help = 'force overwrite',
      default = False, action = 'store_true')
    args=parser.parse_args()
    
    # path normalisation
    import os
    args._in = os.path.realpath(args._in)
    args.ref = os.path.realpath(args.ref)
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheng/%pheno_%maf.%gset.gsa.out',__file__)
    proj.add_output(args._in+'/%pheng/summary/siggenes.csv',__file__)
    try: main(args)
    except: cmdhistory.errlog()