#!/usr/bin/env python3
'''
Creates diagnostic plots for the asymmetry statistics
'''

def main(args):
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import os
    import fnmatch
    import time
    
    os.chdir(args._in)
    files = []
    for f in os.listdir():
      if fnmatch.fnmatch(f,f'*{args.file}*'):
        files.append(f)
    
    if not os.path.isdir(args.out):
      os.system(f'mkdir -p {args.out}/')
    
    tic = time.perf_counter()
    for i in range(len(files)):
      f = files[i]
      print(f'processing {f}')
      if os.path.isfile(f'{args.out}/{f}'.replace('txt',args.fmt)) and (not args.force):
        print(f'{f} is already plotted')
        continue
      df = pd.read_csv(f,sep = ' ')
      cols = df.columns.values
      cols = cols[cols!='FID']
      cols = cols[cols!='IID']
      df = df[cols]
      df = df.melt(value_name = 'value', var_name = 'metric')
      _,ax = plt.subplots(figsize = (15,10))
      sns.histplot(data = df, x = 'value', hue = 'metric',bins = 50, common_bins = False,
                   ax = ax, legend = False, palette = 'pastel6')
      toc = time.perf_counter() - tic
      print(f'output to {args.out}/{f}, time = {toc:.3f}'.replace('txt','png'))
      plt.savefig(f'{args.out}/{f}'.replace('txt',args.fmt))
      toc= time.perf_counter()-tic
      print(f'finished {f}. {i+1}/{len(files)}, {toc:.3f} seconds')
  
  
if __name__ == '__main__':
    import argparse
    # directory is given in the input
    parser = argparse.ArgumentParser(description='This programme creates diagnostic plots for the IDPs')
    parser.add_argument('-i','--in', dest = '_in', default = '../pheno/conn/',
                        help = 'Data directory')
    parser.add_argument('-o','--out',dest  = 'out', default = '../pheno/plots/',
                        help = 'Output directory')
    parser.add_argument('--file',dest = 'file', default = '*.txt', help = 'Files to plot')
    parser.add_argument('--fmt', dest = 'fmt', default = 'png', help = 'Output Format')
    parser.add_argument('-f','--force', dest = 'force', help = 'Force overwrite',
                        default = False, const = True, action = 'store_const')
    args = parser.parse_args()
    import os
    for arg in ['_in','out']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in, __file__)
    proj.add_output(args._out, __file__)
    try: main(args)
    except: cmdhistory.errlog()