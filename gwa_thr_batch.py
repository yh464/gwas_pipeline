#!/usr/bin/env python3

def main(args):
    # modules
    import pandas as pd
    import numpy as np
    import os
    from fnmatch import fnmatch
    
    # parse input
    thr = args.thr.split(',')
    for i in range(len(thr)):
      thr[i] = float(thr[i])
      exec(f'dflist_{i} = []')
    if not os.path.isdir(args.out): os.mkdir(args.out)
    if args.force: f = ' -f'
    else: f = ''
    
    for x in args.pheno:
      os.chdir(args._in)
      os.chdir(x)
      print(f'Processing phenotype {x}')
      if not os.path.isdir(f'{args.out}/{x}/'): os.mkdir(f'{args.out}/{x}/')
      if not os.path.isdir(f'{args.out}/{x}/cache/'): os.mkdir(f'{args.out}/{x}/cache/')
      # initialise lists
      for i in range(len(thr)):
         exec(f'thr{i}_list = []')
      
      # all GWA files
      flist = []
      for y in os.listdir():
        if fnmatch(y, '*.fastGWA') and (not fnmatch(y, '*X.fastGWA')) and \
          (not fnmatch(y, '*all_chrs.fastGWA')):
          flist.append(y)
      
      # generate file-wise siglist
      sigflist = []
      for y in flist:
        prefix = y.replace('.fastGWA','')
        q = True
        tmp = []
        for i in range(len(thr)):
          tmp.append(f'{args.out}/{x}/cache/siglist.{thr[i]}.{prefix}.txt')
          if not os.path.isfile(f'{args.out}/{x}/cache/siglist.{thr[i]}.{prefix}.txt'): q = False
        sigflist.append(tmp)
        if q and (not args.force): continue
        scripts_path = os.path.realpath(__file__)
        scripts_path = os.path.dirname(scripts_path)
        os.system(f'bash {scripts_path}/pymaster.sh '+
          f'gwa_thr.py -i {args._in}/{x}/{y} -t {args.thr} -o {args.out}/{x}/cache/ {f}')
      
      # concatenate siglist (always overwrites)
      sigflist = np.array(sigflist, dtype = 'U')
      for i in range(len(thr)):
        q = True
        for y in sigflist[:,i]:
          if not os.path.isfile(y): q = False
        if not q: print(f'Files not fully generated for {x}, please run again!'); continue
        out_fname = f'{args.out}/{x}/siglist.{thr[i]}.txt'
        with open(out_fname, 'w') as file:
          print('NODE\tCHR\tSNP\tPOS\tA1\tA2\tN\tAF1\tBETA\tSE\tP', file = file)
          file.close()
        
        for y in sigflist[:,i]:
          os.system(f'cat {y} >> {out_fname}')
        
        df = pd.read_csv(out_fname, sep = '\s+')
        df.sort_values(by = ['CHR','POS','P'], inplace = True)
        df.to_csv(out_fname, sep = '\t', index = False)
        exec(f'dflist_{i}.append(df)')
    
    # log output
    for i in range(len(thr)):
      exec(f'tmp = pd.concat(dflist_{i})')
      tmp.sort_values(by = ['CHR','POS','P'], inplace = True)
      tmp.to_csv(f'{args.out}/siglist.{thr[i]}.txt', sep = '\t', index = False)
      tmp = tmp['SNP'].unique()
      with open(f'{args.out}/sigstats.{thr[i]}.txt','w') as file:
        print(f'Number of significant SNPs at threshold {thr[i]}:', file = file)
        print(len(tmp), file = file)
        print('\n\n\n\n\n', file = file)
        for x in tmp: print(x, file = file)
        file.close()
        
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
      description = 'This script filters out the most significant SNPs for regional traits')
    parser.add_argument('pheno', nargs='*', help = 'phenotypes to process',
      default = ['deg_local','degi_local','degc_local','clu_local','eff_local','mpl_local'])
    parser.add_argument('-i','--in', dest = '_in', help = 'GWA master directory',
      default = '../gwa/')
    parser.add_argument('-o','--out',dest = 'out', help = 'output directory')      # defaults to input dir
    parser.add_argument('-t','--thr', dest= 'thr',
      help = 'significance threshold. format: 5e-8,1e-5 (separate different thresholds w/ commas)',
      default = '3.1075e-11,5e-8')
    parser.add_argument('-f','--force', dest = 'force', help = 'force overwrite',
                        action = 'store_true', default = False)
    args = parser.parse_args()
    import os
    args._in = os.path.realpath(args._in)
    if type(args.out) == type(None): args.out = args._in+'/siglist/'
    else: args.out = os.path.realpath(args._out)
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_var('%reg',r'.+', 'brain region')
    proj.add_var('%pheng',r'.+', 'phenotype group')
    proj.add_var('%maf',r'[0-9.]+', 'minor allele freq') # only allows digits and decimals
    proj.add_var('%sig',r'[0-9.+-]+', 'significance level')
    proj.add_input(args._in+'/%pheng/%reg_%maf.fastGWA', __file__)
    proj.add_input(args._in+'/%pheng/%reg_%maf_X.fastGWA', __file__)
    proj.add_output(args.out+'/siglist/%pheng/siglist.%sig.txt', __file__)
    try: main(args)
    except: cmdhistory.errlog()