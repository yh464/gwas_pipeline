#!/usr/bin/env python3
'''
This script concatenates the PRS-cs output to produce individual level PRS scores
'''
    
def main(args):
    import os
    import numpy as np
    import pandas as pd
    
    # parse input
    flist = np.loadtxt(args._list, dtype = 'U')
    bed_list = open(args.bed).read().splitlines()
    if not os.path.isdir(args.out): os.mkdir(args.out)
    
    for i in range(flist.shape[0]):
      f = flist[i,0]
      prefix = f.split('/')[-1].replace('.txt','')
      in_dir = args._in + '/'+ prefix + '/'
      out_dir = args.out +'/'+ prefix +'/'                                     # scores by chr in args.out/prefix/
      if not os.path.isdir(out_dir): os.mkdir(out_dir)
      
      # Score by chromosome
      out_flist = []
      for j in range(22):
        effsz = in_dir + f'chr{j+1}_pst_eff_a1_b0.5_phi{args.phi:.0e}_chr{j+1}.txt'
        out_fname = f'{out_dir}/{prefix}.chr{j+1}'
        out_flist.append(out_fname+'.sscore')
        if not os.path.isfile(out_fname+'.sscore') or args.force:
          os.system(f'{args.plink} --bfile {bed_list[j]} --score {effsz} 2 4 6 center '+
                  f'cols=fid,denom,dosagesum,scoresums --out {out_fname}')
      
      # sum across all chromosomes
      for j in range(22):
        f = out_flist[j]
        df = pd.read_csv(f, sep = '\s+')
        if j == 0: 
          ref = df.iloc[:,:2]
          total_score = df.iloc[:,-1]
        else: 
          df = df.merge(ref, how = 'inner')
          total_score += df.iloc[:,-1]
      
      out = ref
      out.columns = ['FID','IID']
      out['score_total'] = total_score
      out['score_norm'] = total_score/total_score.std()
      out.to_csv(f'{args.out}/{prefix}.txt', index = False, sep = '\t')
     
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'This script concatenates the PRS-cs output to produce individual level PRS scores')
    parser.add_argument('--list', dest = '_list',
      help = 'list of GWA files to generate PRS. FORMAT: ${absolute path} \t ${sample size}',
      default = '../params/gwa_for_prs.list')
    parser.add_argument('-i','--in', dest = '_in', help = 'input directory',
      default = '../prs/prs_effsize/')
    parser.add_argument('--plink', dest = 'plink', help = 'Location of plink2 executable',
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Genetics/plink2')
    parser.add_argument('--bed', dest = 'bed', help = 'list of PLINK binaries for target sample',
      default = '../params/bed_files_ukb.txt')
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
      default = '../prs/prs_score/')
    parser.add_argument('--phi', dest = 'phi', help = 'shrinkage parameter used for prscs',
      type = float, default = 0.01)
    parser.add_argument('--force','-f', dest = 'force', action = 'store_true',
                        default = False, help = 'force overwrite')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','plink','bed','_list']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheno.txt', __file__)
    proj.add_output(args.out+'/%pheno.txt', __file__)
    try: main(args)
    except: cmdhistory.errlog()