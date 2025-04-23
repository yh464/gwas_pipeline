#!/usr/bin/env python3
'''
python script that runs the asym_stats.py script in batches
'''

def main(args):
    import os
    import numpy as np
    import pandas as pd
    import time
    
    # nroi: HCP = 376, 500sym = 334, aparc = 84, economo = 102, sjh = 1027
    nroi = 376
      
    if args.force:
      fs = ' -f'
    else: fs = ''
    
    # array submitter
    from _utils.slurm import array_submitter
    submitter = array_submitter(
        name = 'pheno',
        partition = 'icelake',
        timeout = 5, mode = 'long',
        debug = False
        )
    
    tic = time.perf_counter()
    subjs = np.loadtxt(args.subjs,dtype = 'U')
    nsubj = subjs.size
    logdir = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/logs/'
    logout = 'asym_stats_batch.log'
    f = open(logdir+logout,'w')
    idx = 0
    
    for subj in subjs:
      if idx % 100 == 0:
          toc = time.perf_counter()-tic
          print(f'Subject {idx} / {nsubj}, time = {toc:.3f}')
      idx += 1
      in_fname = args._in.replace('%sub',subj)
      
      # check progress
      if not os.path.isfile(in_fname):
        print(f'{subj}: no connectome found', file = f)
        continue
      
      skip = True
      if args.force: skip = False
      
      if skip:
        try:
          tmp = np.loadtxt(f'{args.out}/global/{subj}.txt')
          if tmp.size != 17: skip = False
        except: skip = False
      
      if skip:
        try:
          tmp = np.loadtxt(f'{args.out}/global_asym/{subj}.txt')
          if tmp.size != 17: skip = False
        except: skip = False
      
      if skip:
        try:
          tmp = pd.read_csv(f'{args.out}/local/{subj}.txt')
          if tmp.shape[0] != nroi or tmp.shape[1] != 7: skip = False
        except: skip = False
      
      if skip:
        try:
          tmp = pd.read_csv(f'{args.out}/local_asym/{subj}.txt')
          if tmp.shape[0] != nroi/2 or tmp.shape[1] != 21: skip = False
        except: skip = False
      
      if skip: 
        print(f'{subj}: connectome is already phenotyped', file = f)
      else:
        indir = args._in.replace('%sub', subj)
        submitter.add(
          f'python pheno.py {subj} -i {indir} -o {args.out}/{fs}')
        print(f'{subj} submitted for analysis', file = f)
    
    submitter.submit()
    
    
if __name__ == '__main__':
    # input argument processing
    import argparse
    from _utils.slurm import parser_config   
    parser = argparse.ArgumentParser(description='This programme processes the connectome '+
                               ' for one single individual for imaging derived phenotypes')
    parser.add_argument('-i','--in',dest = '_in', help =
        'Target file to screen',
        default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Imaging/'+
        'UKB%sub/func/fMRI/parcellations/HCP.fsaverage.aparc_seq/Connectivity_sc2345.txt')
    parser.add_argument('-s', '--subjs', dest = 'subjs', help = 'list of subjs',
        default = '../params/subjlist_rsfmri_hcp.txt')
    parser.add_argument('-o','--out',dest = 'out', help = 'Output directory',
        default = '../pheno/ukb/')
    parser.add_argument('-f','--force', dest = 'force', help = 'Force output',
        default = False,const = True, action = 'store_const')
    
    args = parser.parse_args()
    import os
    for arg in ['_in','out','subjs']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
        
    from _utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in, __file__)
    proj.add_output(args.out, __file__)
    try: main(args)
    except: cmdhistory.errlog()