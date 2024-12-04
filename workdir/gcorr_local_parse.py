#!/usr/bin/env python
'''
This script parses regional correlations and plots heatmaps
'''

def main(args):
    import os
    from fnmatch import fnmatch
    import pandas as pd
    import numpy as np
    from scipy.stats import pearsonr
    from time import perf_counter as t
    
    tic = t()
    
    # parse input
    os.chdir(args._in)
    if args.pheno == 'all':
      tmp = []
      for x in os.listdir():
        if fnmatch(x, '*local'): tmp.append(x)
      args.pheno = tmp
    if type(args.out) == type(None): args.out = args._in
    
    log = open(f'{args.out}/qc.txt','w')
    
    # generate prefix list
    dflist = []
    for x in args.pheno:
      os.chdir(args.sumstats)
      os.chdir(x)
      
      flist = []
      for y in os.listdir():
        if fnmatch(y, '*.sumstats'):
          flist.append(y.replace('.sumstats',''))
      flist = pd.Series(flist, name = 'roi')
      dflist.append(flist)
    
    # prefix is the common subset of all GWA files for each pheno
    flist = dflist[0]
    for x in dflist:
      flist = pd.merge(flist, x)
    flist = flist['roi'].tolist()
    nroi = len(flist)
    nphen = len(args.pheno)
    
    toc = t() - tic
    print(f'Parsed input. time = {toc:.3f}')
    
    # parse global correlations
    def parse_rg(f):
      f = open(f).read().splitlines()
      tmp = False
      for x in f:
        if tmp: dat = x; break                                                     # this takes the next line of h2_obs_se
        if fnmatch(x, '*h2_obs_se*'): tmp = True
      dat = dat.split(' ')
      while True:
        try: dat.remove('')
        except: break
      try: rg = float(dat[2])
      except: rg = np.nan
      try: z = float(dat[4])
      except: z = np.nan
      try: gcov = float(dat[-2])
      except: gcov = np.nan
      return rg, z, gcov
    
    global_rg = np.zeros((nroi,nphen))
    global_z = np.zeros((nroi,nphen))
    global_gcov = np.zeros((nroi,nphen))
    
    i = 0
    for x in args.pheno:
      j = 0
      os.chdir(args._in)
      os.chdir(x)
      os.chdir('idp_rg')
      for y in flist:
        rg,z,gcov = parse_rg(f'global.{y}.rg.log')
        global_rg[j,i] = rg
        global_z[j,i] = z
        global_gcov[j,i] = gcov
        j += 1
      i += 1
    
    global_rg = pd.DataFrame(global_rg,index = flist, columns = args.pheno)
    global_z = pd.DataFrame(global_z,index = flist, columns = args.pheno)
    global_gcov = pd.DataFrame(global_gcov,index = flist, columns = args.pheno)
    global_rg.to_csv(f'{args.out}/global_rg.csv')
    global_z.to_csv(f'{args.out}/global_z.csv')
    global_gcov.to_csv(f'{args.out}/global_gcov.csv')
    
    toc = t() - tic
    print(f'Parsed global genetic correlations. time = {toc:.3f}')
    
    # parse regional correlations
    rglist = []
    zlist = []
    gcovlist = []
    
    for x in args.pheno:
      out_rg = f'{args.out}/{x}_rg.csv'
      out_z = f'{args.out}/{x}_z.csv'
      out_gcov = f'{args.out}/{x}_gcov.csv'
      
      if os.path.isfile(out_rg) and os.path.isfile(out_z) and os.path.isfile(out_gcov) \
        and (not args.force): 
        rglist.append(pd.read_csv(out_rg,index_col = 0))
        zlist.append(pd.read_csv(out_z, index_col = 0))
        gcovlist.append(pd.read_csv(out_gcov, index_col = 0))
        continue
      
      os.chdir(args._in)
      os.chdir(x)  
      os.chdir('idp_rg')
      
      xrg = np.zeros((nroi,nroi))
      xz = np.zeros((nroi,nroi))
      xgcov = np.zeros((nroi,nroi))
      idx = 0
      for i in range(nroi):
        y = flist[i]
        xz[i,i] = 0; xgcov[i,i] = 0                                                # zero diagonals except rg
        for j in range(nroi):
          k = flist[j]
          if i == j: continue
          idx += 1
          try:
            rg,z,gcov = parse_rg(f'{y}.{k}.rg.log')
          except:
            rg,z,gcov = parse_rg(f'{k}.{y}.rg.log')
          xrg[i,j] = rg; xrg[j,i] = rg
          xz[i,j]= z; xz[j,i] = z
          xgcov[i,j]= gcov; xgcov[j,i] = gcov
          
          if idx % 1000 == 0:
            toc = t() - tic
            print(f'Parsed regional genetic correlations ({idx}/{nroi**2}). time = {toc:.3f}')
      
      os.chdir('../h2/')
      for i in range(nroi):
        y = flist[i]
        f = open(f'{y}.h2.log').read().splitlines()
        for k in f:
          if fnmatch(k, 'Total Observed scale h2*'):
            l = k.split(' ')
            break
        xrg[i,i] = float(l[-2])
      
      xrg = pd.DataFrame(xrg,index = flist, columns = flist)
      xgcov = pd.DataFrame(xgcov,index = flist, columns = flist)
      xz = pd.DataFrame(xz,index = flist, columns = flist)
      xrg.to_csv(out_rg)
      xz.to_csv(out_z)
      xgcov.to_csv(out_gcov)
      rglist.append(xrg); zlist.append(xz); gcovlist.append(xgcov)
    
    toc = t() - tic
    print(f'Parsed regional genetic correlations (complete). time = {toc:.3f}')
    
    # quality control
    for x in args.pheno:
      for y in flist:
        if np.isnan(global_rg.loc[y,x]):
          print(f'{args.out}/{x}/idp_rg/global.{y}.rg.log', file = log)
    print('\n'*5,file = log)
    for x in args.pheno:
      for y in flist:
        if np.isnan(global_rg.loc[y,x]):
          print(f'{args.out}/{x}/h2/{y}.h2.log', file = log)
    print('\n'*5,file = log)
    for i in range(nphen):
      for j in range(nroi):
        for k in range(j):
          tmp = rglist[i]
          if np.isnan(tmp.loc[flist[j],flist[k]]):
            if os.path.isfile(f'{args.out}/{args.pheno[i]}/idp_rg/{flist[j]}.{flist[k]}.rg.log'):
              print(f'{args.out}/{args.pheno[i]}/idp_rg/{flist[j]}.{flist[k]}.rg.log', file = log)
            else:
              print(f'{args.out}/{args.pheno[i]}/idp_rg/{flist[k]}.{flist[j]}.rg.log', file = log)
      print('\n'*5, file = log)
      
    # test cophenetic correlation
    def mantelr(m1, m2):
      n = m1.shape[0]
      if m1.shape[1] != n or (m2.shape[0] != n or m2.shape[1] != n):
        raise ValueError('Must compare two SQUARE matrices with the same size')
      lt = np.tri(n, k=-1).astype('?')
      tmp = pd.DataFrame(dict(l = m1.values[lt], r = m2.values[lt])).dropna()
      res = pearsonr(tmp['l'].values,tmp['r'].values)
      return res.statistic
    
    rg_mantel = np.zeros((nphen, nphen))
    for i in range(nphen):
      for j in range(i):
        rg_mantel[i,j] = mantelr(rglist[i],rglist[j])
    rg_mantel += rg_mantel.T
    rg_mantel = pd.DataFrame(rg_mantel, index = args.pheno, columns = args.pheno)
    rg_mantel.to_csv(f'{args.out}/rg_mantel.csv')
    
    # plot figures
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
    
    for i in range(nphen):
      x = args.pheno[i]
      out_hmap = f'{args.out}/{x}_rg.pdf'
      if os.path.isfile(out_hmap) and (not args.force): continue
      _,ax = plt.subplots(figsize = (25,25))
      sns.heatmap(rglist[i],vmin = -1, vmax = 1, cmap = 'redblue',
                  ax = ax)
      plt.savefig(out_hmap, bbox_inches = 'tight')
      plt.close()
    toc = t() - tic
    print(f'Plotted rg figures. time = {toc:.3f}')
    
    for i in range(nphen):
      x = args.pheno[i]
      out_hmap = f'{args.out}/{x}_z.pdf'
      if os.path.isfile(out_hmap) and (not args.force): continue
      _,ax = plt.subplots(figsize = (25,25))
      sns.heatmap(zlist[i],vmin = -1, vmax = 1, cmap = 'redblue',
                  ax = ax)
      plt.savefig(out_hmap, bbox_inches = 'tight')
      plt.close()
    toc = t() - tic
    print(f'Plotted z figures. time = {toc:.3f}')
    
    _,ax = plt.subplots(1,2,figsize = (13,25))
    sns.heatmap(global_rg, vmin = -1, vmax = 1, cmap = 'redblue', ax = ax[0])
    sns.heatmap(global_z, vmin = -2.5, vmax = 2.5, cmap = 'redblue',
                ax = ax[1])
    plt.savefig(f'{args.out}/global_rg_z.pdf', bbox_inches = 'tight')
    plt.close()
    toc = t() - tic
    print(f'Total time {toc:.3f}')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
      description = 'This script parses regional correlations and plots heatmaps')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*', default = 'all')
    parser.add_argument('-i','--in', help = 'input directory', dest = '_in',
      default = '../local_corr/')
    parser.add_argument('-s','--sumstats', help = 'GWA directory', dest = 'sumstats',
      default = '../gene_corr/ldsc_sumstats/') # this is used to parse file names
    parser.add_argument('-o','--out', help = 'output directory', dest = 'out',
      default = None) # defaults to input directory
    parser.add_argument('-f','--force', help = 'force overwrite', action = 'store_true',
      default = False)
    args = parser.parse_args()
    
    import os
    for arg in ['_in','gwa']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    if type(args.out) == type(None): args.out = args._in
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheng/idp_rg/%reg.%reg.rg.log', __file__)
    proj.add_input(args.sumstats+'/%pheng/%pheno_%maf.fastGWA',__file__)
    proj.add_output(args.out+'/%pheno_.*.csv', __file__) # .* is a wildcard
    try: main(args)
    except: cmdhistory.errlog()