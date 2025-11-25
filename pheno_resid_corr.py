#!/usr/bin/env python

def main(args):
    import os
    import pandas as pd
    os.chdir(args._in)
    
    from fnmatch import fnmatch
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
    
    df = pd.read_csv(f'{args.pheno[0]}_resid.txt')
    for i in range(1,len(args.pheno)):
      df = df.merge(pd.read_csv(f'{args.pheno[i]}_resid.txt'), on = ['FID','IID']).dropna()
    
    df = df.iloc[:,2:]
    out_col = []
    for x in df.columns.tolist():
      if (not fnmatch(x, '*frac*')) and (not fnmatch(x, '*abs*') and (not fnmatch(x, 'bet*'))):
        out_col.append(x)
    df = df[out_col]
    
    corr = df.corr()
    out_fname = args.out + '/resid_corr_'+'_'.join(args.pheno)
    corr.to_csv(out_fname + '.csv')
    sns.heatmap(corr, vmin = -1, vmax = 1, cmap = 'redblue', annot = True, fmt = '.2f',
                    ax = plt.subplots(figsize = (22,22))[1])
    
    plt.savefig(out_fname+ '.png')
    plt.close()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Computes correlational plots between IDPs and PGS (single file)')
    parser.add_argument('pheno', nargs = '*', help = 'Phenotypes to correlate',
      default = ['global_norm','global_asym_norm','global_structural'])
    parser.add_argument('-i','--in', dest = '_in', help = 'Input directory',
      default = '../pheno/resid/')
    parser.add_argument('-o','--out',dest = 'out', help = 'Output directory')
    # always forces
    args=parser.parse_args()
    import os
    args._in = os.path.realpath(args._in)
    if type(args.out)==type(None): args.out = args._in
    
    from ._utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in, __file__)
    proj.add_output(args._out, __file__)
    try: main(args)
    except: cmdhistory.errlog()