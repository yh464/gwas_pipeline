#!/usr/bin/env python3
'''
concatenates global, local, and asymmetry IDPs
'''

def main(args):
    import numpy as np
    import pandas as pd
    import os
    import time

    os.chdir(args._in)
    
    # subject list
    os.chdir(f'{args._in}/global/')
    slist = np.array(os.listdir(),dtype = 'U')
    for i in range(slist.size):
      slist[i] = slist[i][3:10]
      
    slist = slist.astype(np.int32)
    slist.sort()
    n = slist.size
    nodes = np.loadtxt(
      '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/params/HCP.fsaverage_4mm_names.txt', 
      dtype ='U')
    nodes_u = nodes[:int(nodes.size/2)].copy()
    for i in range(nodes_u.size):
      nodes_u[i] = nodes_u[i][2:] # unilateral label
    
    # global phenotypes
    tic = time.perf_counter()
    
    if args.force or (not os.path.isfile(f'{args._in}/global.txt')):
      data = np.zeros((17,n))
      for i in range(n):
        s = slist[i]
        tmp = np.loadtxt(f'UKB{s}.txt')
        data[:,i] = tmp.copy()
        if i % 10000 == 0: 
            toc = time.perf_counter()-tic
            print(f'{i}/{n*4}, time = {toc:.3f} seconds')
      for i in range(17):
          data[i,:] /= data[i,:].std() # normalises standard deviation to 1
      df = pd.DataFrame(dict(FID = slist,
                             IID = slist,
                             deg_global = data[0,:],
                             degi_global = data[1,:],
                             degc_global = data[2,:],
                             mpl_global = data[3,:],
                             eff_global = data[4,:],
                             clu_global = data[5,:],
                             smw_global = data[6,:],
                             degi_l = data[7,:],
                             degi_r = data[8,:],
                             mpl_l = data[9,:],
                             mpl_r = data[10,:],
                             eff_l = data[11,:],
                             eff_r = data[12,:],
                             clu_l = data[13,:],
                             clu_r = data[14,:],
                             smw_l = data[15,:],
                             smw_r = data[16,:]))
      os.chdir(args._in)
      df = df.dropna()
      df.to_csv('global.txt',index = False, sep = ' ')
    
    # local phenotypes
    if args.force or (not os.path.isfile(f'{args._in}/bet_local.txt')):
      os.chdir(f'{args._in}/local/')
      data = np.zeros((376,7,n))
      for i in range(n):
        s = slist[i]
        tmp = pd.read_csv(f'UKB{s}.txt',sep = ',')
        data[:,:,i] = tmp.values.copy()
        if i == 0:
          local_columns = tmp.columns
        if i % 10000 == 0: 
            toc = time.perf_counter()-tic
            print(f'{n+i}/{n*4}, time = {toc:.3f} seconds')
      for i in range(376):
          for j in range(7):
              data[i,j,:] /= data[i,j,:].std()
      os.chdir(args._in)
      for i in range(local_columns.size):
        df = pd.DataFrame(data = data[:,i,:].T,columns = nodes)
        df.insert(0,'IID',slist)
        df.insert(0,'FID',slist)
        df.to_csv(f'{local_columns[i]}.txt',index = False, sep = ' ')
        
    # local asymmetry
    if args.force or (not os.path.isfile(f'{args._in}/bet_nasym_rank.txt')):
      os.chdir(f'{args._in}/local_asym/')
      data = np.zeros((188,21,n))
      for i in range(n):
        s = slist[i]
        tmp = pd.read_csv(f'UKB{s}.txt',sep = ',')
        data[:,:,i] = tmp.values.copy()
        if i == 0:
          local_columns = tmp.columns
        if i % 10000 == 0:
            toc = time.perf_counter()-tic
            print(f'{n*2+i}/{n*4}, time = {toc:.3f} seconds')
      for i in range(188):
          for j in range(21):
              data[i,j,:] /= data[i,j,:].std()
      os.chdir(args._in)
      for i in range(local_columns.size):
        df = pd.DataFrame(data = data[:,i,:].T,columns = nodes_u)
        df.insert(0,'IID',slist)
        df.insert(0,'FID',slist)
        df.to_csv(f'{local_columns[i]}.txt',index = False, sep = ' ')
    
    # global asymmetry
    if args.force or (not os.path.isfile(f'{args._in}/global_asym.txt')):
      os.chdir(f'{args._in}/global_asym/')
      data = np.zeros((17,n))
      
      for i in range(n):
        s = slist[i]
        tmp = np.loadtxt(f'UKB{s}.txt')
        data[:,i] = tmp.copy()
        if i % 10000 == 0: 
            toc = time.perf_counter()-tic
            print(f'{n*3+i}/{n*4}, time = {toc:.3f} seconds')
      
      df = pd.DataFrame(dict(FID = slist,
                             IID = slist,
                             degi_asym_abs = data[0,:],
                             degi_asym_frac = data[1,:],
                             mpl_asym_abs = data[2,:],
                             mpl_asym_frac = data[3,:],
                             eff_asym_abs = data[4,:],
                             eff_asym_frac = data[5,:],
                             clu_asym_abs = data[6,:],
                             clu_asym_frac = data[7,:],
                             smw_asym_abs = data[8,:],
                             smw_asym_frac = data[9,:],
                             deg_asym_corr = data[10,:],
                             degi_asym_corr = data[11,:],
                             degc_asym_corr = data[12,:],
                             mpl_asym_corr = data[13,:],
                             eff_asym_corr = data[14,:],
                             clu_asym_corr = data[15,:],
                             con_asym_corr = data[17,:],
                             bet_asym_corr = data[16,:]),index = slist)
      for i in range(17):
          data[i,:] /= data[i,:].std() # normalises standard deviation to 1
      os.chdir(args._in)
      df = df.dropna()
      df.to_csv('global_asym.txt',index = False, sep = ' ')
    
if __name__ == '__main__':
    # input argument processing
    import argparse as ap
    parser = ap.ArgumentParser(description='This programme concatenates all subjects for their functional connectome phenotypes')
    parser.add_argument('-i','--in', dest = '_in', 
                        default = '../pheno/ukb/',
                        help = 'Data directory')
    parser.add_argument('-f','--force', dest = 'force', help = 'Force output',
                        default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    args._in = os.path.realpath(args._in)
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in, __file__)
    proj.add_output(args._in, __file__)
    try: main(args)
    except: cmdhistory.errlog()