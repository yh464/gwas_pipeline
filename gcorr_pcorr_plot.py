#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2023-07-13

Creates a plot of genetic, phenotypic correlations and heritability

Requires following inputs: 
    GWAS summary statistics (scans directory for all files)
'''


def main(args):
    import os
    from fnmatch import fnmatch
    import pandas as pd
    import numpy as np
    
    os.chdir(args._in)
    if type(args.pheno) == type('a'):
      pheno = [args.pheno]                                                         # forcibly convert to list 
    else:
      pheno = args.pheno
    
    # Merge phenotype dataframes
    df_list = []
    pheno_list = []
    for x in pheno:
      df_list.append(pd.read_csv(f'{x}.txt', sep = '\s+'))
      pheno_list += [x] * (len(df_list[-1].columns) - 1)
    pheno_df = df_list[0]
    for x in range(1,len(df_list)):
      pheno_df = pheno_df.merge(df_list[x])
    del df_list
    
    # Generate correlational matrix
    pheno_df = pheno_df.iloc[:,2:]
    n = pheno_df.shape[0]
    prefix_list = pheno_df.columns.tolist()
    tmp = []
    for i in prefix_list:
      if not (fnmatch(i,'*frac*') or fnmatch(i,'*abs*') or i == 'bet_asym_corr'):  # remove absolute and fractional differences
        tmp.append(i)
    prefix_list = tmp
    print(prefix_list)
    
    pheno_df = pheno_df[prefix_list]
    corrmat = pheno_df.corr()
    semat = ((1-corrmat**2)/(n-2)) **0.5                                           # STE of Pearson's rho 
    del pheno_df
    
    # Parse logfiles of heritability
    for x in args.pheno:
      os.chdir(args.h2)
      os.chdir(x)
      ls = os.listdir()
      for y in ls:
        for z in prefix_list:
          if fnmatch(y,f'{z}*.greml.hsq'):
            df = pd.read_csv(y, sep = '\t').iloc[:10,1:].astype(float)
            corrmat.loc[z,z] = df.values[3,0]
            semat.loc[z,z] = df.values[3,1]
            del df
    
    # Parse logfiles of genetic correlation
    os.chdir(args.corr)
    ls = os.listdir()
    for x in range(len(prefix_list)):
      for y in range(x,len(prefix_list)):                                           # upper triangle = genetic
        x1 = prefix_list[x]; x2 = prefix_list[y]
        p1 = pheno_list[x]; p2 = pheno_list[y]
        for z in ls:
          if fnmatch(z,f'{p1}_{x1}*{p2}_{x2}*rg.log') or \
            fnmatch(z,f'{p2}_{x2}*{p1}_{x1}*rg.log'):
            tmp = open(z)
            tmp_stats = tmp.readlines(-1)
            tmp_stats = tmp_stats[-4].split(' ')
            while tmp_stats.count('') > 0: # remove empty fields
              tmp_stats.remove('')
            try: corrmat.iloc[x,y] = float(tmp_stats[2])
            except: 
                corrmat.iloc[x,y] = np.nan
                print(z)
            try: semat.iloc[x,y] = float(tmp_stats[3])
            except: semat.iloc[x,y] = np.nan
    
    # Output csv
    prefix = '_'.join(args.pheno)
    corrmat.to_csv(f'{args.out}/correlation_{prefix}.csv', index = True, header = True)
    semat.to_csv(f'{args.out}/se_{prefix}.csv', index = True, header = True)
    
    # Z score and P values
    zscore = corrmat / semat
    import scipy.stats as sts
    import numpy as np
    pval = sts.norm.cdf(-np.abs(zscore)) * 2
    
    # FDR output
    pfdr = sts.false_discovery_control(pval)
    pfdr = pd.DataFrame(pfdr, columns=semat.columns, index = semat.index)
    pfdr.to_csv(f'{args.out}/pfdr_{prefix}.csv', index = True, header = True)
    
    # Plot heatmap
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
    
    corrmat['id'] = corrmat.index
    summary = pd.melt(corrmat, id_vars = 'id')
    
    # point sizes by significance
    pval = pd.DataFrame(pval, columns = semat.columns, index = semat.index)
    pval_m = pd.melt(pval)
    pfdr_m = pd.melt(pfdr)
    size = np.ones(summary.shape[0])
    for i in range(summary.shape[0]):
        if pval_m.loc[i,'value'] < 0.05: size[i] += 1
        if pfdr_m.loc[i,'value'] < 0.05: size[i] += 1
    summary['size'] = size
    
    # scatter plot style heatmap
    sns.set_theme(style = 'whitegrid')
    _,ax = plt.subplots(figsize = (10.5,10))
    sns.scatterplot(
        summary,
        x = 'id', y = 'variable',
        hue = 'value', palette = 'redblue', hue_norm = (-1,1),
        size = 'size', sizes = (50,200),
        edgecolor = '.7',
        ax = ax, legend = False
        )
    plt.xlabel('phenotypic correlation')
    plt.ylabel('genetic correlation')
    
    # colour bar
    norm = mpl.colors.Normalize(vmin=-1, vmax=1)
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap='redblue'), ax = ax)
    
    # diagonal line
    ax.axline((0,0), slope = 1, color = 'k', zorder = 0)
    ax.annotate('Heritability', (-1,-1),xytext=(-5,-4.5), rotation = 45)
    
    # aesthetics
    for _, spine in ax.spines.items():
      spine.set_visible(False)
    ax.invert_yaxis()
    for label in ax.get_xticklabels():
      label.set_rotation(90)  
    plt.savefig(f'{args.out}/correlation_{prefix}.pdf', bbox_inches = 'tight')
    plt.savefig(f'{args.out}/correlation_{prefix}.png', bbox_inches = 'tight')
    plt.close()
    
    # H2 estimate bar plot
    h2_toplot = np.diag(corrmat)
    se_toplot = np.diag(semat)
    _,ax = plt.subplots(figsize = (12, 6))
    ax.bar(x = prefix_list, height = h2_toplot, yerr = se_toplot)
    ax.tick_params(axis = 'x', labelrotation = 45)
    plt.savefig(f'{args.out}/h2_se_{prefix}.pdf', bbox_inches = 'tight')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'This programme parses genetic correlation outputs and generates heatmaps')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
    parser.add_argument('-i','--in', dest = '_in', help = 'Phenotype file directory',
      default = '../pheno/ukb/')
    parser.add_argument('--corr', dest = 'corr', help = 'correlational logfile directory',
      default = '../gene_corr/gcorr/')
    parser.add_argument('--h2', dest = 'h2', help = 'Heritability logfile directory',
      default = '../gwa/')
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
      default = '../gene_corr/')
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
      default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','corr','h2']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args.corr+'/gcorr/%pheno_%maf.%pheno_%maf.rg.log', __file__)
    proj.add_input(args.h2+'/%pheng/%pheno_%maf.h2.log',__file__)
    proj.add_output(args.out+'/correlation_%pheno..*', __file__) # .* is a wildcard
    proj.add_output(args.out+'/h2_%pheno..*', __file__) # .* is a wildcard)
    try: main(args)
    except: cmdhistory.errlog()