#!/usr/bin/env python3
'''
this programme parses the hsq files from GREML output
'''

def main(args):
    import os
    from fnmatch import fnmatch
    import pandas as pd
    
    if type(args.pheno) == type('a'):
      pheno = [args.pheno]                                                         # forcibly convert to list 
    else:
      pheno = args.pheno
    
    # temp directory
    tmpdir = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/temp'
    if not os.path.isdir(tmpdir): os.system(f'mkdir -p {tmpdir}')
    
    flist = []
    vg = []; ve = []; vp = []; h2 = []; ll = []; ll0 = []; lrt = []; dof = []; p = []; n = []
    se_vg = []; se_ve = []; se_vp = []; se_h2 = []
    
    for x in pheno:
      os.chdir(args._in)
      os.chdir(x)
      for y in os.listdir():
        if fnmatch(y,'*.hsq'):
          df = pd.read_csv(y, sep = '\t').iloc[:10,1:].astype(float)
          vg.append(df.values[0,0]); se_vg.append(df.values[0,1])
          ve.append(df.values[1,0]); se_ve.append(df.values[1,1])
          vp.append(df.values[2,0]); se_vp.append(df.values[2,1])
          h2.append(df.values[3,0]); se_h2.append(df.values[3,1])
          ll.append(df.values[4,0])
          ll0.append(df.values[5,0])
          lrt.append(df.values[6,0])
          dof.append(df.values[7,0])
          p.append(df.values[8,0])
          n.append(df.values[9,0])
          flist.append(y.replace('.greml.hsq',''))
    
    df = pd.DataFrame(dict(
      phenotype = flist,
      h2 = h2,
      h2_se = se_h2,
      p = p,
      vg = vg,
      vg_se = se_vg,
      ve = ve,
      ve_se = se_ve,
      vp = vp,
      vp_se = se_vp,
      ll = ll,
      ll0 = ll0,
      lrt = lrt,
      dof = dof,
      n = n))
    
    df.to_csv(f'{args.out}/greml_summary.csv', index = False)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'This programme parses the hsq files'+
      ' from GREML analysis')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
    parser.add_argument('-i','--in', dest = '_in', help = 'GWA file directory',
      default = '../gwa/')
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
      default = '../gwa/')
    # always overwrites
    args = parser.parse_args()
    import os
    for arg in ['_in','out']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheng/%pheno.greml..*', __file__)
    proj.add_output(args.out+'/greml_summary.csv', __file__)
    try: main(args)
    except: cmdhistory.errlog()