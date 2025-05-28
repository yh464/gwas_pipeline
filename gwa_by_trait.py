#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-02-24

Wrapper for fastGWA analysis for a single phenotype

Requires following inputs: 
    phenotype file in FID IID *** format
    PLINK bed binaries
    covariants files in FID IID *** format
'''

def main(args):
    import pandas as pd
    from fnmatch import fnmatch
    
    # parse bed files
    if fnmatch(args.bed, '*.bed'):
        bfile = f'--bfile {args.bed[:-4]} --autosome'
    elif os.path.isfile(f'{args.bed}.bed'):
        bfile = f'--bfile {args.bed} --autosome'
    elif os.path.isfile(args.bed):
        try: 
            open(args.bed).read().splitlines()
            bfile = f'--mbfile {args.bed} --autosome'
        except:
            raise ValueError('Please supply a valid BED file prefix')
    elif os.path.isdir(args.bed):
        bed_list = []
        temp_blist = f'{args.bed}/bed.list'
        with open(temp_blist,'w') as f:
            for j in range(22):
                for y in os.listdir(args.bed):
                    if fnmatch(y.lower(), f'*chr{j+1}.bed'): 
                        bed_list.append(f'{args.bed}/{y[:-4]}')
                        print(bed_list[-1], file = f)
        bfile = f'--mbfile {temp_blist} --autosome'
    else: raise ValueError('Please supply a valid BED file prefix')
    
    # parse MAF
    if args.maf != 'raw':
      ft = f'--maf {args.maf}'
    else: ft = ''
    
    # read X chromosome bed file
    if fnmatch(args.xbed, '*.bed'):
        xbfile = f'--bfile {args.xbed[:-4]}'
    elif os.path.isfile(f'{args.xbed}.bed'):
        xbfile = f'--bfile {args.xbed}'
    else: args.xchr = False; print('Warning: skipping X chromosome because no bed file found')
    
    if not os.path.isfile(f'{args.out}.fastGWA') or args.force:
        os.system(f'{args.gcta} --fastGWA-mlm {bfile} --grm-sparse {args.grm} '+
            f'--pheno {args._in} --mpheno {args.mpheno} --qcovar {args.qcov} --covar {args.dcov}'+
            f' {ft} --keep {args.keep} --out {args.out}')
    
    if not args.xchr: return
    if not os.path.isfile(f'{args.out}_X.fastGWA') or args.force:
        xkeep_file = args.keep.replace('.txt','_X.txt')
        if not os.path.isfile(xkeep_file):
            xkeep = pd.read_table(args.keep)
            xfam = pd.read_table(f'{xbfile}.fam'.replace('--bfile ',''), header = None, usecols = [0,1])
            xfam.columns = ['FID','IID']
            xkeep = pd.merge(xkeep, xfam)
            xkeep.to_csv(xkeep_file, sep = '\t', index = False)
        
        os.system(f'{args.gcta} --fastGWA-mlm {bfile} --grm-sparse {args.grm} '+
            f'--pheno {args._in} --mpheno {args.mpheno} --qcovar {args.qcov} --covar {args.dcov}'+
            f' {ft} --keep {xkeep_file} --model-only --out {args.out}_Xmodel')
        os.system(f'{args.gcta} {xbfile} --load-model {args.out}_Xmodel.fastGWA --geno 0.1 --out {args.out}_X')
        os.system(f'tail -n +2 {args.out}_X.fastGWA >> {args.out}.fastGWA')

if __name__ == '__main__':
    import argparse
    # argument input
    parser = argparse.ArgumentParser(description=
      'This programme runs GWA for any phenotype given as the 1st positional argument')
    
    io = parser.add_argument_group(title = 'input and output options')
    io.add_argument('-i','--in', dest = '_in', help = 'Phenotype file', required = True)
    io.add_argument('-o','--out',dest  = 'out', help = 'Output prefix', required = True)
    io.add_argument('--mpheno', help = 'Phenotype ID in file, required for fastGWA')
    io.add_argument('--dcov',dest = 'dcov', help = 'DISCRETE covariance file',
      default = '../params/discrete_covars.txt')
    io.add_argument('--qcov',dest = 'qcov', help = 'QUANTITATIVE covariance file',
      default = '../params/quantitative_covars.txt')
    io.add_argument('--bed',dest = 'bed', help = 'PLINK2 binaries',
      default = '../params/bed_files_ukb.txt')
    io.add_argument('--grm', dest = 'grm', help = 'Genetic correlation matrix',
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/sp0.05_grm')
    
    params = parser.add_argument_group(title = 'parameters for GCTA')
    params.add_argument('--gcta', dest = 'gcta', help = 'Location of GCTA executable',
      default = '../toolbox/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1')
    params.add_argument('--maf', dest = 'maf', help = 'Filter by minor allele frequency',
      default = '0.01', type = str)
    params.add_argument('--keep', dest = 'keep', help = 'Subjects to keep', # intentionally absolute
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/ukbkeepfile_202402.txt')
    
    xchr = parser.add_argument_group(title = 'X chromosome GWAS options')
    xchr.add_argument('--nox', dest = 'xchr', help = 'Do not conduct GWAS for X chromosome',
        default = True, action = 'store_false')
    xchr.add_argument('--xbed', help = 'PLINK binary for the X chromosome',
        default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/bed/chrX')
    
    parser.add_argument('-f','--force', dest = 'force', help = 'Force overwrite',
      default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','gcta','dcov','qcov','grm','bed','xbed']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
    from _utils import logger
    logger.splash(args)
    
    main(args)