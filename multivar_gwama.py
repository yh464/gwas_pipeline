#!/usr/bin/env python3
    
def main(args):
    import os
    from fnmatch import fnmatch
    import pickle
    import time
    import pandas as pd
    import numpy as np
    import scipy.stats as sts
    import matplotlib.pyplot as plt
    
    tic = time.perf_counter()
    
    # parse input list
    flist = open(args._list).read().splitlines()
    flist = args.file+flist
    prefix = []
    dirs = []
    
    for x in flist:
      if not fnmatch(x,'*.fastGWA'):
        flist.remove(x)
      tmp = x.split('/')
      tmp = tmp[-1]
      dirs.append(x.replace(tmp,''))                                               # remove the directories
      tmp = tmp.replace('.fastGWA', '')
      prefix.append(tmp)
    n = len(flist)
    
    # output directory
    if not os.path.isdir(args.out): os.mkdir(args.out)
    os.chdir(args.out)
    
    # log and temp
    log = open(f'/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/logs/gwama_{args.prefix}.log','w')
    print('Traits analysed in n-weighted GWAMA', file = log)
    for x in prefix:
      print(x, file = log)
    tmp = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/temp/'
    if not os.path.isdir(tmp): os.mkdir(tmp)
    
    # extract data
    print(file = log)
    print('Extracting data for analysis', file = log)
    
    # load cache
    cache = f'{tmp}gwama_{args.prefix}_cache.pickle'
    if os.path.isfile(cache) and (not args.force):
      for k, v in pickle.load(open(cache,'rb')).items():
        globals()[k] = v # dflist, ref, h2, cti
      toc = time.perf_counter()-tic
      print(f'Loaded prepared data from cache. time = {toc:.3f} seconds', file = log)
      
    # compile cache
    else:
      dflist = []
      cti = np.diag(np.ones(n)) # cross-trait intercept
      h2 = np.zeros(n)
      
      for i in range(n):
        x = flist[i]
        prefix_x = prefix[i]
        
        # read fastGWA file
        df = pd.read_csv(x, sep = '\s+').drop_duplicates(subset = 'SNP')
        df.insert(loc = df.shape[1]-1, column = 'Z', value = df.BETA/df.SE)
        # columns: chr, snp, pos, a1, a2, n, af1, beta, se, z, p
        dflist.append(df[['N','BETA','SE','Z','P']])
        
        # read h2 file
        if not os.path.isfile(x.replace('_0.01.fastGWA','.greml.hsq')): raise ValueError
        f = open(x.replace('_0.01.fastGWA','.greml.hsq')).read().splitlines()
        for y in f:
          if 'V(G)/Vp' in y:
            h2[i] = float(y.split('\t')[1])
            break
        
        # read rg file
        for j in range(i):
          prefix_y = prefix[j]
          try:
            f = open(f'{args.rg}/{prefix_x}.{prefix_y}.rg.log').read().splitlines()
          except:
            f = open(f'{args.rg}/{prefix_y}.{prefix_x}.rg.log').read().splitlines()
          for k in range(-1, -len(f),-1):
            if 'gcov_int' in f[k]:
              v = f[k+1].split(' ') # gcov_int is in the next line
              break
          while True: # remove all blank instances
            try: v.remove('')
            except: break
          if float(v[-2]) > 1: v[-2] = 1
          if float(v[-2]) < -1: v[-2] = -1
          cti[i,j] = float(v[-2])
          cti[j,i] = float(v[-2])
        
        toc = time.perf_counter() - tic
        print(f'Processed {prefix_x} ({i+1}/{n}), time = {toc:.3f} seconds')
        print(f'Processed {prefix_x} ({i+1}/{n}), time = {toc:.3f} seconds', file = log)
      
      ref = df[['CHR', 'SNP', 'POS', 'A1', 'A2', 'AF1']]
      
      # write cache
      with open(cache, 'wb') as f:
        pickle.dump(dict(
          dflist = dflist,
          ref = ref,
          h2 = h2,
          cti = cti), f)
    
    print(h2)
    print(cti)
    # We skip the sanity checks since all GWA data have been generated in the same pipeline
    
    # Extract weights and weighted Z-score by trait
    w = []
    nsnp = ref.shape[0]
    wz_total = np.zeros(nsnp)
    neff = np.zeros(nsnp)
    print('Calculating Neff for each SNP x trait', file = log)
    print('Calculating Neff for each SNP x trait')
    for i in range(n):
      neff += dflist[i]['N']
      weight = dflist[i]['N'] **0.5 * h2[i]**0.5
      w.append(weight.tolist())
      dflist[i]['wZ'] = weight * dflist[i]['Z']
      print(dflist[i].head())
      wz_total += dflist[i]['wZ']
    
    print('Calculating the weighted z-score for each SNP', file = log)
    print('Calculating the weighted z-score for each SNP')
    w = np.array(w) # n rows, nsnp columns
    print(w.shape)
    for i in range(nsnp):
      t = w[:,i].reshape((n,1))
      coef = (t.T @ cti @ t)
      if coef == 0:
        print(f'WARNING: {ref.SNP[i]} has zero coefficient')
      wz_total[i] /= coef**0.5
    
    print('Calculating the p-values for each SNP', file = log)
    print('Calculating the p-values for each SNP')
    p_total = 1-sts.chi2.cdf(wz_total**2, df = 1)
    
    # We skip the AF1 calculation because the sample is totally overlapping
    
    # calculate beta and SE where se sqrt(1/neff^2 / maf / (1-maf))
    print('Calculating beta values for each SNP')
    print('Calculating beta values for each SNP', file = log)
    beta_total = wz_total / neff / (ref.AF1 * (1-ref.AF1))**0.5
    se_total = beta_total / wz_total
    out = ref
    out.insert(5, column = 'N', value = neff)
    out['BETA'] = beta_total
    out['SE'] = se_total
    out['P'] = p_total
    
    # write output
    toc = time.perf_counter() - tic
    print(f'Writing file to {args.out}/{args.prefix}.gwama, time = {toc:.3f} seconds', file = log)
    out.to_csv(f'{args.out}/{args.prefix}.gwama', sep = '\t', index = False)
    
    toc = time.perf_counter() - tic
    msg = f'Analysis finished at {toc:.3f} seconds. {nsnp} SNPs were included'
    print(msg, file = log)
    print(msg)
    
    # Manhattan and qqplot
    from qmplot import manhattanplot, qqplot
    if (not os.path.isfile(f'{args.out}/{args.prefix}.gwama.png')) or args.force:
      manhattanplot(data = df,
                  chrom = 'CHR',
                  pos = 'POS',
                  pv = 'P',
                  snp = 'SNP',
                  logp = True,
                  text_kws = {'fontfamily','monospace'})
      plt.savefig(f'{args.out}/{args.prefix}.gwama.png')
      plt.close()
    
    if (not os.path.isfile(f'{args.out}/{args.prefix}.gwama.qqplot.png')) or args.force:
      qqplot(data =df['P'], title = x.replace('.fastGWA',''),
             marker= '.', xlabel=r"Expected $-log_{10}{(P)}$",
               ylabel=r"Observed $-log_{10}{(P)}$")
      plt.savefig(f'{args.out}/{args.prefix}.gwama.qqplot.png')
      plt.close()
    toc = time.perf_counter()-tic
    print(f'Fig plotted, time = {toc:.3f} seconds.', file = log)
    print(f'Fig plotted, time = {toc:.3f} seconds.')
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'This script constitutes the '
      +'GWAMA pipeline for a given list of GWA files')
    parser.add_argument('--list', dest = '_list', help = 'List of fastGWA files to include')
    parser.add_argument('--file', dest = 'file', help = 'Any additional fastGWA files',
                        nargs = '*', default = [])
    parser.add_argument('--rg', dest = 'rg', help = 'Directory of LDSC rg log files',
      default = '../gene_corr/gcorr/')
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
      default = '../multivar-gwa/')
    parser.add_argument('-p','--prefix', dest = 'prefix', help = 'name of the output file',
                        required = True)
    parser.add_argument('-f', '--force', dest = 'force', action = 'store_true',
                        default = False, help = 'force overwrite')
    args = parser.parse_args()
    import os
    for arg in ['_list','out','rg']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    
    from ._utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._list, __file__)
    proj.add_output(args.out+'/%pheng.gwama.fastGWA', __file__)
    try: main(args)
    except: cmdhistory.errlog()