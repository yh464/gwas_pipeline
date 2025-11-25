def main(args):
    import os
    import pandas as pd
    from fnmatch import fnmatch
    
    if not os.path.isdir(args.out): os.mkdir(args.out)
    
    # load input files
    os.chdir(args._in)
    flist = []
    for g in args.pheno:
        for f in os.listdir(f'{args._in}/{g}'):
            if fnmatch(f, '*.txt'): flist.append(f'{args._in}/{g}/{f}')
    
    # merge datasets
    df = pd.read_table(args.ref).drop('FID', axis = 'columns')
    for f in flist:
        prefix = f[:-4] # remove .txt
        tmp = pd.read_table(f).drop(['FID','score_norm'], axis = 'columns')
        tmp.columns = ['IID',prefix]
        df = pd.merge(df, tmp)
    
    # correlation matrix
    df.index = df.IID
    df = df.drop('IID', axis = 'columns')
    corr = df.corr()
    corr.to_csv(f'{args.out}/prs_corr_qc.txt', sep = '\t', header = True, index = True)
    
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
  
    _,ax = plt.subplots(figsize = (5,5))
    sns.heatmap(corr, vmin = -1, vmax = 1, cmap = 'redblue', ax = ax)
    plt.savefig(f'{args.out}/prs_corr_qc.pdf', bbox_inches = 'tight')
    plt.savefig(f'{args.out}/prs_corr_qc.png', bbox_inches = 'tight')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
      description = 'this script quality controls PRS scores')
    parser.add_argument('pheno', help = 'Phenotype groups to generate PRS',
      nargs = '*', default = ['disorders','disorders_subtypes'])
    parser.add_argument('-i','--in', dest = '_in', help = 'input directory',
      default='../prs/prs_score/')
    parser.add_argument('-r','--ref',dest = 'ref', help = 'reference standard UKB PRS',
      default = '/home/yh464/rds/rds-rb643-ukbiobank2/Data_Users/yh464/params/ukb_std_prs.txt')
    # intended to be absolute path
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory')

    # always overwrites
    args = parser.parse_args()
    import os
    args._in = os.path.realpath(args._in)
    if type(args.out) == type(None): args.out = f'{args._in}/qc/'
    
    from ._utils import cmdhistory
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()