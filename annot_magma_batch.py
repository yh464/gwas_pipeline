#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2024-11-12

Conducts MAGMA for single GWAS summary statistics
To be used in combination with annot_magma_by_trait.py

Requires following inputs: 
    GWAS summary statistics (single file)
    MAGMA annotation files (including H-MAGMA)
    MAGMA binary
'''

def main(args):
    import os
    from fnmatch import fnmatch
    
    # array submitter
    from _utils.slurm import array_submitter
    submitter = array_submitter(
        name = 'annot',
        timeout = 240,
        debug = False
        )
    
    # parse input
    f = open(args.gset)
    glist = f.readlines()
    if args.force: force = '-f'
    else: force = ''
    logdir = os.path.realpath('../logs/')
    if not os.path.isdir(logdir): os.mkdir(logdir)
    if not os.path.isdir(args.out): os.mkdir(args.out)
    
    for x in args.pheno:
      os.chdir(args._in)
      os.chdir(x)
      if not os.path.isdir(f'{args.out}/{x}'): os.mkdir(f'{args.out}/{x}')
      
      # find all files with only autosomes
      flist = []
      for y in os.listdir():
        if fnmatch(y,'*.fastGWA') and (not fnmatch(y,'*_X.fastGWA')) and (not fnmatch(y, '*_all_chrs*')):
          flist.append(y)
      
      for y in flist:
        skip = True
        prefix = y.replace('.fastGWA','')
        out_geneset = f'{args.out}/{x}/{prefix}'
        
        # check existing files
        if not os.path.isfile(f'{out_geneset}.genes.raw'): skip = False
        for z in glist:
          z = z.replace('\n', '')
          zprefix = z.split('/')[-1].replace('.txt','')
          gsaout = f'{args.out}/{x}/{prefix}.{zprefix}'
          if (not os.path.isfile(f'{gsaout}.gsa.out')): skip = False
        
        if skip and (not args.force): continue
        
        # batch run
        scripts_path = os.path.realpath(__file__)
        scripts_path = os.path.dirname(scripts_path)
        submitter.add('python '+
          f'annot_magma_by_trait.py {x} -i {y} -d {args._in} -a {args.annot} -b {args.bfile} '+
          f'-o {args.out} --magma {args.magma} --gset {args.gset} {force}')
    submitter.submit()

if __name__ == '__main__':
    import argparse
    from _utils.slurm import parser_config
    parser = argparse.ArgumentParser(
      description = 'This file batch runs the MAGMA pipeline for all IDPs')
    parser.add_argument('pheno', nargs = '*', help = 'Phenotypes')
    parser.add_argument('-i','--in', dest = '_in', help = 'Directory containing all phenotypes',
      default = '../gwa/')
    parser.add_argument('-a', '--annot', dest ='annot', help = 'directory to annotation files',
      default = '../toolbox/hmagma')
    parser.add_argument('-b', '--bfile', dest = 'bfile', help = 'bed binary to use in magma',
      default = '../toolbox/magma/g1000_eur')
    parser.add_argument('-o', '--out', dest = 'out', help = 'output directory',
      default = '../annot/magma')
    parser.add_argument('--magma', dest = 'magma', help = 'MAGMA executable',
      default = '../toolbox/magma/magma')
    parser.add_argument('--gset', dest ='gset', help = 'Gene sets to study enrichment',
      default = '../params/gset.txt')
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
      default = False, action = 'store_true')
    parser = parser_config(parser)
    args = parser.parse_args()
    
    # path normalisation
    import os
    for arg in ['_in','annot','bfile','out','magma','gset']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheng/%pheno_%maf.fastGWA', __file__)
    proj.add_output(args.out+'/%pheng/%pheno_%maf.log',__file__)
    try: main(args)
    except: cmdhistory.errlog()