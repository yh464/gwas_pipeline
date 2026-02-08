
import os
import pandas as pd
import numpy as np
import scipy.stats as sts
from tqdm import tqdm
from _utils.path import find_gwas
from _utils.plugins.logparser import crosscorr_parse
from _utils import logger
from multiprocessing import Pool, cpu_count
log = logger.logger()

def read_sumstats(input_args):
    g, p, in_dir, weight = input_args
    df = pd.read_table(f'{in_dir}/{g}/{p}.fastGWA', usecols = 
        lambda x: x.upper() in (['CHR','SNP','POS','A1','A2','BETA','OR','SE','N', 'AF1']),
        index_col = ['CHR','SNP','POS','A1','A2'])
    df.columns = [x.upper() for x in df.columns]
    if 'OR' in df.columns and not 'BETA' in df.columns: df['BETA'] = np.log(df['OR'])
    
    n = (df['N']).rename((g,p))
    w = (weight * (n ** 0.5)).rename((g,p))
    z = (df['BETA'] * w / df['SE']).rename((g,p))
    af = (df['AF1'] * n).rename((g,p))
    return w, z, af, n

def main(args):
    pheno = find_gwas(args.pheno, dirname = args._in, long = True)

    log.log('Estimating the weights for each trait')
    corr = crosscorr_parse(pheno, full = True)
    rg = corr.pivot(index = ['group1','pheno1'], columns = ['group2','pheno2'], values = 'rg')
    gcovint = corr.pivot(index = ['group1','pheno1'], columns = ['group2','pheno2'], values = 'gcov_int')
    if args.pca: # use PCA to estimate weights
        pc1 = np.linalg.eig(rg.fillna(0).values)[1][:,0]
        weights = pd.Series(pc1, index = rg.index)
    elif args.nw: # use n-weighted meta-analysis
        weights = pd.Series(np.diag(rg.fillna(0).values), index = rg.index) ** 0.5
    else:
        raise ValueError('Please specify a method to estimate weights: --pca or --nw')
    
    log.log('This script assumes all alleles are in the same order across all files. Please run gwa_harmonise.py before calling this script.')

    # weight_list = []
    # aflist = []
    # idx = pd.MultiIndex.from_frame(pd.DataFrame(index = [], columns = ['CHR','SNP','POS','A1','A2']))
    # n_total = pd.Series(name = 'N', index = idx, dtype = float)
    # out_z = pd.Series(name = 'Z', index = idx, dtype = float)
    # for g, p in tqdm(pheno, desc = 'Reading summary statistics and aggregating weighted Z-scores'):
    #     df = pd.read_table(f'{args._in}/{g}/{p}.fastGWA', usecols = 
    #         lambda x: x.upper() in (['CHR','SNP','POS','A1','A2','BETA','OR','SE','N', 'AF1']),
    #         index_col = ['CHR','SNP','POS','A1','A2'])
    #     df.columns = [x.upper() for x in df.columns]
    #     if 'OR' in df.columns and not 'BETA' in df.columns: df['BETA'] = np.log(df['OR'])

    #     n_total = n_total.add(df['N'], fill_value = 0)
    #     weight = weights.loc[(g,p)] * (df['N'] ** 0.5).rename((g,p))
    #     out_z = out_z.add(df['BETA'] * weight / df['SE'], fill_value = 0) # weighted z-score
    #     weight_list.append(weight)
    #     aflist.append(df['AF1'].rename((g,p)) * df['N'])
    
    # parallelised version of the above loop
    parallel_args = [(g, p, args._in, weights.loc[(g,p)]) for g, p in pheno]
    with Pool(cpu_count() * 4) as pool:
        out = list(tqdm(pool.imap(read_sumstats, parallel_args), total = len(parallel_args), desc = 'Reading summary statistics and aggregating weighted Z-scores'))
    weight_list, z_list, af_list, n_list = zip(*out)
    n_total = pd.concat(n_list, axis = 1).fillna(0).sum(axis = 1)
    out_z = pd.concat(z_list, axis = 1).fillna(0).sum(axis = 1)

    # adjust AF1
    log.log('Estimating allele frequencies weighted by sample size')
    af_list = pd.concat(af_list, axis = 1).fillna(0)
    af1 = af_list.sum(axis = 1).rename('AF1') / n_total

    # for each SNP, divide the weighted Z-score by sqrt(weight[:,SNP].T dot gcov_int dot weight[:,SNP])
    weight = pd.concat(weight_list, axis = 1).fillna(0)
    gcovint = gcovint.loc[weight.columns, weight.columns].fillna(0).values
    div_coef = ((weight.values @ gcovint) * weight.values).sum(axis = 1) ** 0.5
    out_z = out_z / div_coef
    out = pd.concat([af1, out_z, n_total], axis = 1)
    out['P'] = sts.norm.sf(abs(out_z)) * 2
    out['BETA'] = out_z / n_total / (af1 * (1-af1)) ** 0.5
    out['SE'] = out['BETA'] / out_z
    out.to_csv(args.out, sep = '\t', index = True, header = True)

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description = 'A wrapper script to conduct genomic PCA or vanilla nw-GWAMA meta-analysis')
    parser.add_argument('pheno', nargs = '+', help = 'Phenotype groups to include in the analysis')
    parser.add_argument('-i','--in', dest = '_in', default = '../gwa', help = 'Directory of input GWAS summary statistics')
    parser.add_argument('--pca', action = 'store_true', help = 'Use PCA to estimate weights')
    parser.add_argument('--nw', action = 'store_true', help = 'Use n-weighted meta-analysis to estimate weights')
    parser.add_argument('-o','--out', help = 'Output file name', required = True)
    args = parser.parse_args()
    
    for arg in ['_in', 'out']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    if not args.pca and not args.nw:
        args.pca = True # default to PCA if no method specified
        log.warn('No method specified for estimating weights, defaulting to PCA (--pca)')
        
    from _utils import cmdhistory
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()