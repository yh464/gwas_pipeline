#!/usr/bin/env python3
'''
Merges and munges alleles for the metaanalytic data
'''

def main(args):
    import os
    from fnmatch import fnmatch
    
    # array submitter
    from _utils import array_submitter
    submitter = array_submitter.array_submitter(
        name = 'gwama_munge',
        timeout = 10, lim = 1)
    scripts_path = os.path.realpath(__file__)
    scripts_path = os.path.dirname(scripts_path)
    
    for x in args.pheno:
      os.chdir(args._in)
      os.chdir(x)
      flist = []
      for y in os.listdir():
        if fnmatch(y, '*.txt'):
          flist.append(y)
      
      for y in flist:
        prefix = y.replace('.txt', '')
        if os.path.isfile(f'{args.out}/{prefix}.sumstats') and (not args.force):
          continue
        submitter.add(f'bash {scripts_path}/pymaster.sh '+
          '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/'+
          f'scripts/ldsc_master.sh munge_sumstats.py --sumstats {y} --out {args.out}/{prefix}'\
          + f' --merge-alleles {args.snp} --chunksize 500000' \
          )
    submitter.submit()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = \
      'This programme correlates the imaging phenotypes w/ extant GWAS metaanalyses sumstats')
    parser.add_argument('pheno', nargs = '*', help = 'Phenotypes of metaanalyses') # DO NOT USE THIS FOR DISORDERS
    parser.add_argument('-i','--in', dest = '_in', help = 'Input directory',
      default = '../gwa/')
    parser.add_argument('--snp', dest = 'snp', help = 'SNPs used for IDP GWAS',
      default = '../params/full.snplist')
    parser.add_argument('-o','--out', dest = 'out', help = 'Output directory',
      default = '../gene_corr/ext_sumstats/')
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
      default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','snp']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheng/%pheno.fastGWA', __file__)
    proj.add_output(args.out+'/%pheno.sumstats', __file__)
    try: main(args)
    except: cmdhistory.errlog()