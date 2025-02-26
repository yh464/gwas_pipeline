#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2023-07-20
Version 2: 2024-11-14

Summarises gene-set level enrichment for HMAGMA and MAGMA outputs

Preceding workflow:
    annot_batch.py
Requires following inputs:
    MAGMA GSA outputs
Changelog:
    Added the 'annot' paramter to reflect the multitude of Hi-C and nearest-gene-based annotations
    Changed the heatmap to a scatterplot-style heatmap
'''

def main(args):
    import os
    import time
    tic = time.perf_counter()
    from fnmatch import fnmatch
    import pandas as pd
    import scipy.stats as sts
    import numpy as np
    
    toc = time.perf_counter()-tic
    print(f'Loaded modules. Time = {toc:.3f} seconds')
    
    # for plotting
    import seaborn as sns
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from _plots import corr_heatmap
    from _utils.path import normaliser
    norm = normaliser()
    
    # annotation methods
    annot_list = []
    for x in os.listdir(args.annot):
        if fnmatch(x, '*.genes.annot'): 
            annot_list.append(x)
    
    # cell types of interest
    f = open(args.gset).read().splitlines()
    gset = []
    for x in f:
      gset.append(x.split('/')[-1].replace('.txt',''))
    
    idx = 0
    overall_fig = []
    for x in args.pheno:
      idx += 1
      os.chdir(args._in)
      os.chdir(x)
      if not os.path.isdir('summary'): os.mkdir('summary')
      
      # analyse for each annotation separately
      for annot in annot_list:
          # prefix for phenotypes
          pflist = []
          for y in os.listdir('genes'):
            if fnmatch(y, '*.genes.raw'):
                y = y.replace('.genes.raw','')
                if fnmatch(y, f'*_{annot}'):
                    y = y.replace(f'_{annot}','')
                    pflist.append(y)

          for prefix in pflist:
            # Merge cell types
            out_fname = f'summary/{prefix}.{annot}.gsasummary.txt'
            dflist = []
            for xprefix in gset:
              tmp = pd.read_csv(f'gsa/{prefix}.{annot}.{xprefix}.gsa.out', sep = '\s+', comment = '#')
              tmp.insert(0,'gene_set',value = xprefix)
              tmp.insert(0,'annot', value = annot)
              tmp.insert(0,'pheno', value = prefix)
              
              # FDR correction
              tmp1 = tmp.P.values.copy()    
              fdr = sts.false_discovery_control(tmp1)
              tmp['Pfdr']= fdr
              
              dflist.append(tmp)
            df = pd.concat(dflist, axis = 'index')
            df.to_csv(out_fname, sep = '\t', index = False)
        
      toc = time.perf_counter()-tic
      print(f'FDR correction completed ({idx}/{len(args.pheno)}). Time = {toc:.3f} seconds')
        
      os.chdir('summary')
      all_annot = []
      for annot in annot_list:
          dflist = []
          for y in sorted(os.listdir()):
            if not fnmatch(y, f'*{annot}.gsasummary.txt'): continue
            tmp = pd.read_table(y)
            dflist.append(tmp)
          summary = pd.concat(dflist)
          all_annot.append(summary)
          
          if not fnmatch(annot,'*ENSG*'): continue # no need to plot for anything other than ENSG
          
          # plotting manipulations may affect the tabular output
          summary = summary.copy()
          summary['group'] = x
          summary['pheno'] = summary.pheno.str.replace('_0.01','')
          summary = summary.rename(columns = {'VARIABLE':'cell type','pheno':'phenotype','P':'p','Pfdr':'q'})
          
          # plot heatmap and miami-like figure for each annotation
          summary['pt_size'] = (summary.q < 0.05).astype(float) + \
              (summary.p < 0.05).astype(float) + 1
          gs = summary.gene_set.unique()
          wratio = [len(summary.loc[summary.gene_set==z, 'cell type'].unique()) for z in gs]
          summary['-log(fdr)'] = -np.log10(summary.q) * (summary['BETA'] > 0)
          fig1, ax1 = plt.subplots(1, len(gs), width_ratios = wratio,
              figsize = (len(summary['cell type'].unique()), 3))
          
          for i in range(len(gs)):
              tmp = summary.loc[summary.gene_set == gs[i],:]
              sns.barplot(tmp, x = 'cell type', y = '-log(fdr)', hue = 'phenotype', ax = ax1[i], legend = False)
              if i > 0: ax1[i].set_ylabel('')
              ax1[i].axhline(-np.log10(0.05), color = 'k')
              ax1[i].axhline(np.log10(0.05), color = 'k')

          fig1.savefig(f'../{annot}.enrichment.barplot.pdf', bbox_inches = 'tight')
          plt.close()
          
          fig = corr_heatmap(summary[['group','phenotype','gene_set','cell type','BETA','p','q']])
          fig.savefig(f'../{annot}.enrichment.pdf', bbox_inches = 'tight')
          plt.close()
          
          overall_fig.append(summary[['group','phenotype','gene_set','cell type','BETA','p','q']])
          
      all_annot = pd.concat(all_annot)
      norm.normalise(all_annot).to_csv(f'{args._in}/{x}/all_enrichment_summary.txt', sep = '\t', index = False)
      toc = time.perf_counter()-tic
      print(f'FINISHED {idx}/{len(args.pheno)}. Time = {toc:.3f} seconds')
     
    fig = corr_heatmap(pd.concat(overall_fig))
    fig.savefig(f'{args._in}/ENSG_enrichment_'+'_'.join(args.pheno)+'.pdf', bbox_inches = 'tight')
      
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'Parses MAGMA GSA outputs for a group of phenotypes')
    parser.add_argument('pheno', nargs = '*', help = 'Phenotypes')
    parser.add_argument('-i','--in', dest = '_in', help = 'MAGMA output directory',
      default = '../annot/magma')
    parser.add_argument('-a', '--annot', dest ='annot', help = 'directory to annotation files',
      default = '../toolbox/hmagma')
    parser.add_argument('--gset', dest ='gset', help = 'cell types to study enrichment',
      default = '../params/gset.txt')
    # no need for 'force'
    args = parser.parse_args()
    # path normalisation
    args.pheno.sort()
    import os
    for arg in ['_in','gset']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_var('%gset',r'.+','gene set')
    proj.add_input(args._in+'/%pheng/genes/%pheno_%maf.genes.outM ', __file__)
    proj.add_output(args._in+'/%pheng/gsa/%pheno_%maf.%gset.gsa.out',__file__)
    proj.add_output(args._in+'/%pheng/%pheno_%maf.gsasummary.txt',__file__)
    try: main(args)
    except: cmdhistory.errlog()