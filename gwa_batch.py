#!/usr/bin/env python3
'''
runs GWAS on traits
'''

def main(args):
    # locate phenotype file
    import os
    import fnmatch
    os.chdir(args._in)
    flist = []
    for f in os.listdir():
      if fnmatch.fnmatch(f,f'*{args.pheno}*') and not(os.path.isdir(f)):                  # search for all files matching args.pheno
        flist.append(f)
    if len(flist) != 1: raise ValueError('Please give only ONE phenotype file')
    
    # array submitter
    from _utils import array_submitter
    submitter = array_submitter.array_submitter(
        name = f'gwa_{args.pheno}',
        timeout = 90, debug = True
        )
    
    # temp and log
    tmpdir = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/temp/'                                 # temporatory dir
    if not os.path.isdir(tmpdir): os.mkdir(tmpdir)
    logout = open('/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/logs/gwa_by_trait.log','w')
    
    # check validity of the phenotype file
    import pandas as pd
    f = flist[0]
    if not fnmatch.fnmatch(f,'*.txt'): raise ValueError('Phenotype file should be in TXT format')
    os.system(f'head {f} -n 5 > {tmpdir}/temp.txt')
    df = pd.read_csv(f'{tmpdir}/temp.txt',sep = '\s+')
    c = df.columns.values
    if c[0] != 'FID' or c[1] != 'IID':
      raise ValueError('Phenotype file should be in the format: FID IID *pheno')
    
    # create output folder
    outdir = f'{args.out}/{f}_{args.maf}/'.replace('.txt','')
    print(outdir)
    if not os.path.isdir(outdir):
      os.system(f'mkdir -p {outdir}')                                              # this also generates args.out
    
    # phenotypes to be analysed
    c = c[2:]
    print('Following traits are to be GWA-analysed:', file = logout)
    for i in c: print(i, file = logout)
    
    # for each phenotype
    for i in range(c.size):
      mpheno = i+1
      trait = c[i]
      out_fname = outdir + trait + f'_{args.maf}'
      skip = False
      # check existing files
      if os.path.isfile(f'{out_fname}.fastGWA'):
        skip = True
        print(f'Trait already analysed for: {trait}', file = logout)
      
      if skip and (not args.force):
        continue
      
      if args.maf != 'raw':
        ft = f'--maf {args.maf}'
      else: ft = ''
      
      scripts_path = os.path.realpath(__file__)
      scripts_path = os.path.dirname(scripts_path)
      submitter.add('bash '+
        f'{scripts_path}/gwa_by_trait.sh {args.gcta}'+
        f' {args.mb} {args.grm} {args._in}/{f} {mpheno} {args.qcov} {args.dcov}'+
        f' {ft} {out_fname}')
    
    submitter.submit()

if __name__ == '__main__':
    import argparse
    # argument input
    parser = argparse.ArgumentParser(description=
      'This programme runs GWA for any phenotype given as the 1st positional argument')
    parser.add_argument('pheno', help = 'Phenotype file in TXT format - please supply ONLY ONE')
    parser.add_argument('-i','--in', dest = '_in', help = 'Phenotype directory',
      default = '../pheno/ukb/')
    parser.add_argument('--gcta', dest = 'gcta', help = 'Location of GCTA executable',
      default = '../toolbox/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1')
    parser.add_argument('-o','--out',dest  = 'out', help = 'Output directory',
      default = '../gwa/')
    parser.add_argument('--maf', dest = 'maf', help = 'Filter by minor allele frequency',
      default = 'raw', type = str)
    parser.add_argument('--dcov',dest = 'dcov', help = 'DISCRETE covariance file',
      default = '../params/discrete_covars.txt')
    parser.add_argument('--qcov',dest = 'qcov', help = 'QUANTITATIVE covariance file',
      default = '../params/quantitative_covars.txt')
    parser.add_argument('--grm', dest = 'grm', help = 'Genetic correlation matrix',
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/'+
      #'Data_Genetics/Genetic_data/Neuroimaging_samples/sp_grm')
      'Data_Users/yh464/params/sp0.05_grm')
    parser.add_argument('--mb',dest = 'mb', help = 'List of PLINK2 files',
      default = '../params/bed_files_ukb.txt')
    parser.add_argument('-f','--force', dest = 'force', help = 'Force overwrite',
      default = False, const = True, action = 'store_const')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','gcta','dcov','qcov','grm','mb']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
    from _utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_var('/%pheng',r'.+', 'phenotype group')
    proj.add_var('/%pheno',r'.+', 'phenotype')
    proj.add_var('/%maf',r'[0-9.]+', 'minor allele freq') # only allows digits and decimals
    proj.add_input(args._in+'/%pheng.txt', __file__)
    proj.add_output(args.out+'/%pheng/%pheno_%maf.fastGWA', __file__)
    try: main(args)
    except: cmdhistory.errlog()