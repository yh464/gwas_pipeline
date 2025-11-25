#!/usr/bin/env python3
'''
This script filters all GWA files for a phenotype
'''

def main(args):
    import os
    from fnmatch import fnmatch
    
    if type(args.pheno) == type('a'):
      pheno = [args.pheno] # forcibly convert to list 
    else:
      pheno = args.pheno
    
    # array submitter
    from ._utils.slurm import array_submitter
    submitter = array_submitter(
        name = 'gwa_filter',
        timeout = 10, partition = 'sapphire')
    
    for x in pheno:
      out_dir = f'{args._in}/{x}_{args.freq}/'
      out_dir = f'{args._in}/{x}/'
      print(f'Output to {out_dir}')
      if not os.path.isdir(out_dir): os.system(f'mkdir -p {out_dir}')
      os.chdir(args._in)
      os.chdir(x+'_raw')
      for y in os.listdir():
        if fnmatch(y,'*.fastGWA'):
          out_fname = y.replace('_raw','')
          submitter.add(f'python gwa_filter.py -i {args._in}/{x}_raw/{y} -o {out_dir}/{out_fname} --freq {args.freq}')
    submitter.submit()
    
if __name__ == '__main__':
    from ._utils.slurm import slurm_parser

    # argument input
    parser = slurm_parser(description=
      'This programme filters the GRM to different thresholds')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
    parser.add_argument('--freq', dest = 'freq', help = 'Minor allele frequency threshold',
                        default = 0.01)
    parser.add_argument('-i','--in', dest = '_in', help = 'GWA file directory',
      default = '../gwa/')
    parser.add_argument('-f','--force', dest = 'force', help = 'Force output',
      default = False, const = True, action = 'store_const')
    
    args = parser.parse_args()
    
    import os
    args._in = os.path.realpath(args._in)
    
    from ._utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_var('%pheng',r'.+', 'phenotype group')
    proj.add_var('%pheno',r'.+', 'phenotype')
    proj.add_var('%maf',r'[0-9.]+', 'minor allele freq') # only allows digits and decimals
    proj.add_input(args._in+'/%pheng_raw/%pheno_raw.fastGWA', __file__)
    proj.add_output(args._in+'/%pheng/%pheno.fastGWA', __file__)
    try: main(args)
    except: cmdhistory.errlog()