#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Computes correlational plots between IDPs and UKB phenotypes
'''

def main(args):
    import os
    from fnmatch import fnmatch
    # from sklearn.linear_model import LinearRegression
    
    # array submitter
    from _utils import array_submitter
    submitter = array_submitter.array_submitter(
        name = 'pheno',
        timeout = 5, lim = 100)
    scripts_path = os.path.realpath(__file__)
    scripts_path = os.path.dirname(scripts_path)
    
    # input parsing
    f = ' -f' if args.force else ''
    if not os.path.isdir(args.out): os.mkdir(args.out)
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
      submitter.add(f'bash {scripts_path}/pymaster.sh '+
        f'ukb_corr_plot.py -i {args._in}/{x} -u {args.ukb} --cov {args.cov} --qcov {args.qcov} -o {args.out} {f}')
    submitter.submit()
      
if __name__ == '__main__':  
    import argparse
    parser = argparse.ArgumentParser(description='Computes corre lational plots between IDPs and PGS (batch)')
    parser.add_argument('-i','--in', dest = '_in', help = 'Input file or directory',
      default = '../pheno/conn/')
    parser.add_argument('-u','--ukb', dest = 'ukb', help = 'UKB phenotype directory',
      default = '../pheno/ukb/')
    parser.add_argument('--cov',dest = 'cov', help = 'DISCRETE covariance file',
      default = '../params/discrete_covars.txt')
    parser.add_argument('--qcov',dest = 'qcov', help = 'QUANTITATIVE covariance file',
      default = '../params/quantitative_covars.txt')
    parser.add_argument('-o','--out', dest = 'out', help = 'Output directory',
      default = '../gene_corr/pcorr/')
    parser.add_argument('-f','--force', dest = 'force', help = 'force overwrite',
                        action = 'store_true', default = False)
    args=parser.parse_args()
    import os
    for arg in ['_in','out','ukb','cov','qcov']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
        
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheno.txt', __file__)
    proj.add_input(args.ukb, __file__)
    proj.add_output(args._out, __file__)
    try: main(args)
    except: cmdhistory.errlog()