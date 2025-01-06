#!/usr/bin/env python3
'''
Computes correlational plots between IDPs and PGS
'''

def main(args):
    import os
    from fnmatch import fnmatch
    # from sklearn.linear_model import LinearRegression
    
    # array submitter
    from _utils import array_submitter
    submitter = array_submitter.array_submitter(
        name = 'prs_corr_plot', n_cpu = 1,
        timeout = 30, lim = 1)
    scripts_path = os.path.realpath(__file__)
    scripts_path = os.path.dirname(scripts_path)
    
    # input parsing
    f = ' -f' if args.force else ''
    if not os.path.isdir(args.out): os.mkdir(args.out)
    diagdir = args.out + '/diagnostics/'
    if not os.path.isdir(diagdir): os.mkdir(diagdir)
    if os.path.isfile(args._in): ls = [args._in]
    elif os.path.isdir(args._in):
      os.chdir(args._in)
      ls = []
      for x in os.listdir():
        if fnmatch(x, '*.txt'):
         if not fnmatch(x, '*nasym*'):
          ls.append(x)
    print('Files to process:')
    for x in ls: print(x)
    
    # Correlation
    for x in ls:
      # os.system(f'sbatch -p cclake -N 1 -n 1 -c 1 -t 00:30:00 -o {logdir}prs.{x}.log '+
      submitter.add(f'bash {scripts_path}/pymaster.sh '+
        f'prs_corr_plot.py -i {args._in}/{x} -p {args.prs} --cov {args.cov} --qcov {args.qcov} -o {args.out} {f}')
    submitter.submit()
    # submitter.debug()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Computes correlational plots between IDPs and PGS (batch)')
    parser.add_argument('-i','--in', dest = '_in', help = 'Input file or directory',
      default = '../pheno/ukb/')
    parser.add_argument('-p','--prs', dest = 'prs', help = 'PRS score directory',
      default = '../prs/prs_score/')
    parser.add_argument('--cov',dest = 'cov', help = 'DISCRETE covariance file',
      default = '../params/discrete_covars.txt')
    parser.add_argument('--qcov',dest = 'qcov', help = 'QUANTITATIVE covariance file',
      default = '../params/quantitative_covars.txt')
    parser.add_argument('-o','--out', dest = 'out', help = 'Output directory',
      default = '../prs/prs_corr/')
    parser.add_argument('-f','--force', dest = 'force', help = 'force overwrite',
                        action = 'store_true', default = False)
    args=parser.parse_args()
    import os
    for arg in ['_in','out','prs','cov','qcov']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
        
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheng.txt', __file__)
    proj.add_output(args.out+'/%pheno_.*.csv', __file__)
    try: main(args)
    except: cmdhistory.errlog()