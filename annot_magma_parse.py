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
    # Red-blue colour map
    cdict = dict(red = ((0,0,0),(1/2,1,1),(1,.8,.8)),
                 green = ((0,0,0),(1/2,1,1),(1,0,0)),
                 blue = ((0,.8,.8),(1/2,1,1),(1,0,0)))
    cmap_name = 'redblue'
    cmap = mpl.colors.LinearSegmentedColormap(cmap_name,cdict,1024)
    try:
      mpl.colormaps.register(cmap)
    except:
      pass
    sns.set_theme(style = 'whitegrid')
    
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
              dflist.append(tmp)
            df = pd.concat(dflist, axis = 'index')
            # FDR correction
            tmp = df.P.values.copy()    
            fdr = sts.false_discovery_control(tmp)
            df['Pfdr']= fdr
            df.to_csv(out_fname, sep = '\t', index = False)
        
      toc = time.perf_counter()-tic
      print(f'FDR correction completed ({idx}/{len(args.pheno)}). Time = {toc:.3f} seconds')
        
      os.chdir('summary')
      all_annot = []
      for annot in annot_list:
          dflist = []
          for y in os.listdir():
            if not fnmatch(y, f'*{annot}.gsasummary.txt'): continue
            tmp = pd.read_table(y)
            dflist.append(tmp)
          summary = pd.concat(dflist)
          all_annot.append(summary)
          
          if not fnmatch(annot,'*ENSG*'): continue # no need to plot for anything other than ENSG
          
          summary['pheno'] = summary.pheno.str.replace('_0.01','')
          summary = summary.rename(columns = {'VARIABLE':'cell type','pheno':'phenotype'})
          
          # plot heatmap and miami-like figure for each annotation
          summary['pt_size'] = (summary.Pfdr < 0.05).astype(float) + \
              (summary.P < 0.05).astype(float) + 1
          beta_max = max([abs(summary['BETA'].max()), abs(summary['BETA'].min())])
          gs = summary.gene_set.unique()
          wratio = [len(summary.loc[summary.gene_set==z, 'cell type'].unique()) for z in gs]
          fig, ax = plt.subplots(1,len(gs),
              width_ratios = wratio,
              figsize = (len(summary['cell type'].unique())/2,len(summary.phenotype.unique())/2))
          
          summary['-log(fdr)'] = -np.log10(summary.Pfdr) * (summary['BETA'] > 0)
          fig1, ax1 = plt.subplots(1, len(gs), width_ratios = wratio,
              figsize = (len(summary['cell type'].unique()), 3))
          
          for i in range(len(gs)):
              tmp = summary.loc[summary.gene_set == gs[i],:]
              sns.scatterplot(
                  tmp,
                  y = 'phenotype', x = 'cell type',
                  hue = 'BETA', palette = 'redblue', hue_norm = (-beta_max, beta_max),
                  size = 'pt_size', sizes = (50, 400),
                  edgecolor = '.7',
                  legend = False, ax = ax[i]
                  )
              for _, spine in ax[i].spines.items():
                    spine.set_visible(False)
              for label in ax[i].get_xticklabels():
                    label.set_rotation(90)
              ax[i].set_xlabel(gs[i])
              if i > 0:
                  ax[i].set_yticklabels([''] * len(ax[i].get_yticklabels()))
              ax[i].set_ylabel('')
              
              sns.barplot(tmp, x = 'cell type', y = '-log(fdr)', hue = 'phenotype', ax = ax1[i], legend = False)
              if i > 0: ax1[i].set_ylabel('')
              ax1[i].axhline(-np.log10(0.05), color = 'k')
              ax1[i].axhline(np.log10(0.05), color = 'k')

          norm = mpl.colors.Normalize(vmin=-beta_max, vmax=beta_max)
          fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap='redblue'), cax = fig.add_axes((0.92, 0.25, 0.02, 0.50)))
          fig.savefig(f'../{annot}.enrichment.pdf', bbox_inches = 'tight')
          fig1.savefig(f'../{annot}.enrichment.barplot.pdf', bbox_inches = 'tight')
          plt.close()  
          
      all_annot = pd.concat(all_annot)
      all_annot.to_csv(f'{args._in}/{x}/all_enrichment_summary.txt', sep = '\t', index = False)
      toc = time.perf_counter()-tic
      print(f'FINISHED {idx}/{len(args.pheno)}. Time = {toc:.3f} seconds')
      
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