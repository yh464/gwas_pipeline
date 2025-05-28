#!/usr/bin/env python3
'''
python script that runs the asym_stats.py script in batches
'''

def main(args):
  import numpy as np
    
  if args.force:
    fs = ' -f'
  else: fs = ''
  
  # array submitter
  from _utils.slurm import array_submitter
  submitter = array_submitter(
    name = 'pheno',
    partition = 'icelake',
    timeout = 5, mode = 'long',
    )
  
  subjs = np.loadtxt(args.subjs,dtype = 'U')
  for subj in subjs:
    subj = subj.replace('UKB','')
    in_fname = args._in.replace('%sub',subj)
    submitter.add(f'python pheno.py {subj} -i {in_fname} -o {args.out} {fs}')
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