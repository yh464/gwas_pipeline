#!/usr/bin/env python3
'''
This programme merges the X chromosome and autosomes
'''

def main(args):
    import os
    from fnmatch import fnmatch
    if len(args.pheno) == 0:
        args.pheno = os.listdir(args._in)
    
    for x in args.pheno:
      os.chdir(args._in)
      if not os.path.isdir(x): continue # sanity check especially for default y = ls(args._in)
      os.chdir(x)
      
      alist = []
      for y in os.listdir():
        if fnmatch(y, '*.fastGWA') and (not fnmatch(y, '*_X.fastGWA')) and \
          (not fnmatch(y, '*_all_chrs.fastGWA')):
          alist.append(y)
      
      for y in alist:
        out_fname = y.replace('.fastGWA', '_all_chrs.fastGWA')
        x_fname = y.replace('.fastGWA','_X.fastGWA')
        os.system(f'cat {y} > {out_fname}')
        os.system(f'tail -n +2 {x_fname} >> {out_fname}')
        
if __name__ == '__main__':
    import argparse
    # argument input
    parser = argparse.ArgumentParser(description=
      'This programme merges the X chromosome and autosomes')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
    parser.add_argument('--freq', dest = 'freq', help = 'Minor allele frequency threshold')
    parser.add_argument('-d','--dir', dest = '_in', help = 'GWA file directory',
      default = '../gwa/')
    # parser.add_argument('-f','--force', dest = 'force', help = 'Force output',
    #   default = False, const = True, action = 'store_const')
    # This programme is too brief that it always forces
    args = parser.parse_args()
    
    import os
    args._in = os.path.realpath(args._in)
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in, __file__)
    proj.add_output(args.out, __file__)
    try: main(args)
    except: cmdhistory.errlog()