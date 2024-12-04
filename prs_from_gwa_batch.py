#!/usr/bin/env python3
'''
This script generates PRScs from external sumstats
'''

def main(args):
    import os
    import numpy as np
    
    # array submitter
    from _utils import array_submitter
    submitter = array_submitter.array_submitter(
        name = 'prs_from_gwa', n_cpu = 2,
        timeout = 150, lim = 1)
    scripts_path = os.path.realpath(__file__)
    scripts_path = os.path.dirname(scripts_path)
    
    # parse input
    flist = np.loadtxt(args._list, dtype = 'U')
    bed_list = open(args.bed).read().splitlines()
    if not os.path.isdir(args.out): os.mkdir(args.out)
    
    for i in range(flist.shape[0]):
      f = flist[i,0]
      n = flist[i,1]
      prefix = f.split('/')[-1].replace('.txt','')
      out_dir = args.out+'/'+prefix+'/'
      if not os.path.isdir(out_dir): os.mkdir(out_dir)
      
      for j in range(22):
        out_prefix = out_dir + f'chr{j+1}'
        if os.path.isfile(out_prefix+f'_pst_eff_a1_b0.5_phi{args.phi:.0e}_chr{j+1}.txt') and (not args.force):
          continue
        submitter.add(f'bash {scripts_path}/prs_from_gwa.sh '+
          f'{args.prscs} {args.ref} {bed_list[j]} {f} {n} {out_prefix} {j+1} {args.phi}')
    submitter.submit()
        
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'This script generates PRS by continuous shrinkage from external sumstats')
    parser.add_argument('--list', dest = '_list',
      help = 'list of GWA files to generate PRS. FORMAT: ${absolute path} \t ${sample size}',
      default = '../params/gwa_for_prs.list')
    parser.add_argument('--prscs', dest = 'prscs', help = 'directory of PRSCS executable',
      default = '/rds/user/yh464/hpc-work/conda/PRScs/')
    parser.add_argument('--ref', dest = 'ref', help = 'reference panel',
      default = '../params/ldblk_1kg_eur/')
    parser.add_argument('--bed', dest = 'bed', help = 'list of PLINK binaries for target sample',
      default = '../params/bed_files_ukb.txt')
    parser.add_argument('-o','--out', dest ='out', help = 'output directory',
      default = '../prs/prs_effsize/')
    parser.add_argument('--phi', dest = 'phi', 
      help = 'shrinkage parameter, set to 10**-4 to 10**-2',
      type = float, default = 0.01)
    parser.add_argument('--force','-f', dest = 'force', action = 'store_true',
                        default = False, help = 'force overwrite')
    args = parser.parse_args()
    import os
    for arg in ['_list','out','prscs','bed','ref']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
        
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._list, __file__)
    proj.add_output(args.out+'/%pheno/.*', __file__)
    try: main(args)
    except: cmdhistory.errlog()