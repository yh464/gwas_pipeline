#!/usr/bin/env python3
'''
Screens faulty files
'''

def main(args):
    import numpy as np
    import os
    import time
    
    subjs = np.loadtxt(args.subjs,dtype = 'U')
    if len(subjs[0]) > 10:
      for i in range(subjs.size):
        subjs[i] = subjs[i][13:23]
    naflag = []
    logout = open('/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/logs/asym_stats_timeout.log','w')
    
    tic = time.perf_counter()
    idx = 0
    for subj in subjs:
      toc = time.perf_counter()-tic
      print(f'{idx}/{subjs.size}, time = {toc:.3f}')
      idx += 1
      # specify directories
      in_filename = args._in.replace('%subj',subj)
      
      if not os.path.isfile(in_filename):
        print(f'No connectome found for the subject {subj}', file = logout)
        continue
      
      # loading the connectome  
      tic = time.perf_counter()
      rsc = np.loadtxt(in_filename,delimiter = ',')
      if np.isnan(rsc).any():
        naflag.append(True)
      else:
        naflag.append(False)

      
    with open(args.out,'w') as f:
      nalist = subjs[naflag]
      for i in nalist:
        print (i, file = f)
    
if __name__ == '__main__':
    # input argument processing
    import argparse as ap
    parser = ap.ArgumentParser(description='This programme processes the connectome '+
                               ' for one single individual for imaging derived phenotypes')
    parser.add_argument('-i','--in',dest = '_in', help =
        'Target file to screen',
        default = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Imaging/'+
        '%subj/func/fMRI/parcellations/HCP.fsaverage.aparc_seq/Connectivity_sc2345.txt')
    parser.add_argument('-s', '--subjs', dest = 'subjs', help = 'list of subjs',
                        default = '../params/subjlist_rsfmri_hcp.txt')
    parser.add_argument('-o','--out',dest = 'out', help = 'Output directory',
                        default = '../params/subjlist_rsfmri_hcp_nan.txt')
    parser.add_argument('-f','--force', dest = 'force', help = 'Force overwrite',
                        default = False,const = True, action = 'store_const')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','subjs']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    
    from ._utils import cmdhistory
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()