#!/usr/bin/env python3
'''
batch runs GREML to assess heritability
'''

def main(args):
    # locate phenotype file
    import os
    import fnmatch
    os.chdir(args._in)
    flist = []
    for f in os.listdir():
      if fnmatch.fnmatch(f,f'*{args.pheno}*') and not(os.path.isdir(f)): # search for all files matching args.pheno
        flist.append(f)
    if len(flist) != 1: raise ValueError('Please give only ONE phenotype file')
    
    # array submitter
    from _utils import array_submitter
    submitter = array_submitter.array_submitter(
        name = f'greml_{args.pheno}', n_cpu = 32,
        timeout = 360, lim = 1)
    
    # temp and log
    tmpdir = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/temp/' # temporatory dir
    if not os.path.isdir(tmpdir): os.mkdir(tmpdir)
    logout = open('/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/logs/gwa_by_trait.log','w')
    
    # check validity of the phenotype file
    import pandas as pd
    f = flist[0]
    if not fnmatch.fnmatch(f,'*.txt'): raise ValueError('Phenotype file should be in TXT format')
    os.system(f'head {f} -n 5 > {tmpdir}/temp.txt')
    df = pd.read_csv(f'{tmpdir}/temp.txt',sep = ' ')
    c = df.columns.values
    if c[0] != 'FID' or c[1] != 'IID':
      raise ValueError('Phenotype file should be in the format: FID IID *pheno')
    
    # create output folder
    outdir = f'{args.out}/{f}/'.replace('.txt','')
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
      out_fname = outdir + trait
      skip = False
      # check existing files
      if os.path.isfile(f'{out_fname}.greml.hsq'):
        skip = True
        print(f'Trait already analysed for: {trait}', file = logout)
      
      if skip and (not args.force):
        continue
     
      scripts_path = os.path.realpath(__file__)
      scripts_path = os.path.dirname(scripts_path)
      submitter.add('bash '+
        f'{scripts_path}/greml_by_trait.sh {args.gcta}'+
        f' {args.grm} {args._in}/{f} {mpheno} {args.cov} {args.qcov}'+
        f' {out_fname}.greml --threads 5')
    
    submitter.submit()
    # submitter.debug()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description=
      'This programme runs GREML for any phenotype given as the 1st positional argument')
    parser.add_argument('pheno', help = 'Phenotype file in TXT format - please supply ONLY ONE')
    parser.add_argument('-i','--in', dest = '_in', help = 'Phenotype directory',
      default = '../pheno/conn/')
    parser.add_argument('--gcta', dest = 'gcta', help = 'Location of GCTA executable',
      default = '../toolbox/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1')
    parser.add_argument('-o','--out',dest  = 'out', help = 'Output directory',
      default = '../gwa/')
    parser.add_argument('--cov',dest = 'cov', help = 'DISCRETE covariate file',
      default = '../params/discrete_covars.txt')
    parser.add_argument('--qcov',dest = 'qcov', help = 'QUANTITATIVE covariate file',
      default = '../params/quantitative_covars.txt')
    parser.add_argument('--grm', dest = 'grm', help = 'Genetic correlation matrix',
      # default = '/rds/project/rb643-1/rds-rb643-ukbiobank2/'+
      # 'Data_Genetics/Genetic_data/Neuroimaging_samples/full_grm')
      # 'Data_Users/yh464/params/sp0.05_grm')
      default = '/rds/project/rds-Nl99R8pHODQ/UKB/Imaging_genetics/yg330/GRM_chr_merged/full_grm')
    parser.add_argument('--mb',dest = 'mb', help = 'List of PLINK2 files',
      default = '../params/bed_files_ukb.txt')
    parser.add_argument('-f','--force', dest = 'force', help = 'Force output',
      default = False, const = True, action = 'store_const')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','gcta','cov','qcov','grm','mb']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheng.txt', __file__)
    proj.add_output(args.out+'/%pheng/%pheno.greml..*', __file__)
    try: main(args)
    except: cmdhistory.errlog()