#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2026-02-08

Python implementation of genomic PCA (Fürtjes et al 2023 https://doi.org/10.1002/hbm.26283) and 
vanilla n-weighted GWAMA (Baselmans et al 2019 https://doi.org/10.1038/s41588-018-0320-8)

Requires following inputs: 
    Harmonised GWAS summary statistics in fastGWA format (CHR, SNP, POS, A1, A2, BETA/OR, SE, N, AF1)
    genetic correlation and heritability estimates (gcorr_batch.py output)
'''

import os, gc
import pandas as pd
import numpy as np
import scipy.stats as sts
from tqdm import tqdm
from _utils.path import find_gwas
from _utils.plugins.logparser import crosscorr_parse
from _utils import logger
from multiprocessing import Pool
log = logger.logger()

def read_snpinfo(input_args):
    g, p, in_dir = input_args
    df = pd.read_table(f'{in_dir}/{g}/{p}.fastGWA', usecols = 
        lambda x: x.upper() in (['CHR','SNP','POS','A1','A2']),
        index_col = ['SNP'], dtype = {
            'CHR': 'category', 'POS': np.int32, 'SNP': str, 'A1': str, 'A2': str
        })
    df = df.loc[~df.index.duplicated(keep = False), :].sort_index()
    df.columns = [x.upper() for x in df.columns]
    return df

def read_sumstats(input_args):
    g, p, in_dir, weight, use_snp = input_args
    df = pd.read_table(f'{in_dir}/{g}/{p}.fastGWA', usecols = 
        lambda x: x.upper() in (['SNP','BETA','OR','SE','N','AF1']),
        index_col = ['SNP'], dtype = {
            'SNP': str, 'BETA': np.float32, 'OR': np.float32, 'SE': np.float64, 'N': np.float32, 'AF1': np.float32
        })
    df = df.loc[~df.index.duplicated(keep = False), :].sort_index()
    df = df.loc[df.index.intersection(use_snp),:]
    df.columns = [x.upper() for x in df.columns]
    if 'OR' in df.columns and not 'BETA' in df.columns: df['BETA'] = np.log(df['OR'])
    n = df['N']
    w = (weight * (n ** 0.5))
    z = (df['BETA'] * w / df['SE'])
    af = (df['AF1'] * n)
    return w, z, af, n

def read_sumstats_snpinfo(input_args):
    g, p, in_dir, weight, = input_args
    df = pd.read_table(f'{in_dir}/{g}/{p}.fastGWA', usecols = 
        lambda x: x.upper() in (['CHR','SNP','POS','A1','A2','BETA','OR','SE','N','AF1']),
        index_col = ['SNP'], dtype = {
            'CHR': 'category', 'POS': np.int32, 'SNP': str, 'A1': str, 'A2': str,
            'BETA': np.float32, 'OR': np.float32, 'SE': np.float64, 'N': np.float32, 'AF1': np.float32
        })
    df = df.loc[~df.index.duplicated(keep = False), :].sort_index()
    df.columns = [x.upper() for x in df.columns]
    if 'OR' in df.columns and not 'BETA' in df.columns: df['BETA'] = np.log(df['OR'])
    n = df['N']
    w = (weight * (n ** 0.5))
    z = (df['BETA'] * w / df['SE'])
    af = (df['AF1'] * n)
    return w, z, af, n, df[['CHR','POS','A1','A2']]

@log.profile
def main(args):
    pheno = find_gwas(args.pheno, dirname = args._in, long = True)

    log.log('Estimating the weights for each trait')
    corr = crosscorr_parse(pheno, full = True)
    rg = corr.pivot(index = ['group1','pheno1'], columns = ['group2','pheno2'], values = 'rg')
    h2 = corr.loc[(corr['group1'] == corr['group2']) & (corr['pheno1'] == corr['pheno2']), ['group1','pheno1','rg']].set_index(['group1','pheno1'])['rg'].rename('h2')
    pheno = h2.loc[h2 > 0].index.tolist() # only include traits with positive heritability
    log.log('Excluding following phenotypes because of negative heritability estimates:')
    for g, p in h2.loc[h2 <= 0].index.tolist():
        log.log(f'    {g} {p} (h2 = {h2.loc[(g,p)]:.4f})')
    rg = rg.loc[pheno, pheno]
    gcovint = corr.pivot(index = ['group1','pheno1'], columns = ['group2','pheno2'], values = 'gcov_int').loc[pheno, pheno]
    if args.pca: # use PCA to estimate weights
        pc1 = np.real(np.linalg.eig(rg.fillna(0).values)[1][:,0])
        weights = pd.Series(pc1, index = rg.index)
    elif args.nw: # use n-weighted meta-analysis
        weights = pd.Series(np.diag(rg.fillna(0).values), index = rg.index) ** 0.5
    else:
        raise ValueError('Please specify a method to estimate weights: --pca or --nw')
    
    log.log('This script assumes all alleles are in the same order across all files. Please run gwa_harmonise.py before calling this script.')

    # read SNP info for each trait
    if args.low_memory:
        parallel_args = [(g, p, args._in) for g, p in pheno]
        with Pool(min(16, len(parallel_args))) as pool:
            snpinfo_list = list(tqdm(pool.imap(read_snpinfo, parallel_args), 
                total = len(parallel_args), 
                desc = 'Reading variant information from summary statistics'))
        
            log.log('Merging variant information across all summary statistics')
            snpinfo = pd.concat(snpinfo_list, axis = 0).drop_duplicates().sort_index()
            if snpinfo.duplicated(['CHR','POS']).any():
                log.warn(f'{snpinfo.duplicated(["CHR","POS"]).sum()} variants have different alleles across files')
                snpinfo.loc[snpinfo.duplicated(['CHR','POS'], keep = False), :].to_csv(args.out.replace('.fastGWA', '.missnp'), sep = '\t', index = False)
                log.log(f'List of variants with inconsistent alleles across files written to {args.out.replace(".fastGWA", ".missnp")}')
                snpinfo = snpinfo.drop_duplicates(['CHR','POS'], keep = False)
                log.log(f'{snpinfo.shape[0]} variants with consistent alleles across files retained for meta-analysis')
            del snpinfo_list
            gc.collect()

            # read summary stats and aggregate weighted Z-scores
            log.log('Reading summary statistics and aggregating weighted Z-scores')
            parallel_args = [(g, p, args._in, weights.loc[(g,p)], snpinfo.index) for g, p in pheno]
            out = list(tqdm(pool.imap(read_sumstats, parallel_args), 
                total = len(parallel_args), 
                desc = 'Reading summary statistics and aggregating weighted Z-scores'))
        weight_list, z_list, af_list, n_list = zip(*out)
        del out

    else:
        parallel_args = [(g, p, args._in, weights.loc[(g,p)],) for g, p in pheno]
        with Pool(min(16, len(parallel_args))) as pool:
            out = list(tqdm(pool.imap(read_sumstats_snpinfo, parallel_args), 
                total = len(parallel_args), 
                desc = 'Reading summary statistics and aggregating weighted Z-scores'))
        weight_list, z_list, af_list, n_list, snpinfo_list = zip(*out)
        log.log('Merging variant information across all summary statistics')
        snpinfo = pd.concat(snpinfo_list, axis = 0).drop_duplicates().sort_index()
        if snpinfo.duplicated(['CHR','POS']).any():
            log.warn(f'{snpinfo.duplicated(["CHR","POS"]).sum()} variants have different alleles across files')
            snpinfo.loc[snpinfo.duplicated(['CHR','POS'], keep = False), :].to_csv(args.out.replace('.fastGWA', '.missnp'), sep = '\t', index = False)
            log.log(f'List of variants with inconsistent alleles across files written to {args.out.replace(".fastGWA", ".missnp")}')
            snpinfo = snpinfo.drop_duplicates(['CHR','POS'], keep = False)
            log.log(f'{snpinfo.shape[0]} variants with consistent alleles across files retained for meta-analysis')
        del snpinfo_list
        gc.collect()

    log.log('Calculating meta-analytic Z-scores')
    n_total = pd.concat(list(n_list), axis = 1, ignore_index = True).fillna(0).sum(axis = 1).rename('N')
    del n_list
    out_z = pd.concat(list(z_list), axis = 1, ignore_index = True).fillna(0).sum(axis = 1).rename('Z')
    del z_list

    # adjust AF1
    log.log('Estimating allele frequencies weighted by sample size')
    af1 = (pd.concat(list(af_list), axis = 1, ignore_index = True).fillna(0).sum(axis = 1) / n_total).rename('AF1')
    del af_list

    # for each SNP, divide the weighted Z-score by sqrt(weight[:,SNP].T dot gcov_int dot weight[:,SNP])
    weight = pd.concat(list(weight_list), axis = 1, ignore_index = True).fillna(0)
    weight.columns = pd.MultiIndex.from_tuples(pheno, names = ['group','pheno'])
    del weight_list
    gcovint = gcovint.loc[weight.columns, weight.columns].fillna(0).values
    div_coef = np.einsum('ij, jk, ik -> i', weight.values, gcovint, weight.values, optimize = 'optimal') ** 0.5
    out_z = out_z / div_coef
    out = pd.concat([snpinfo, af1, out_z, n_total], axis = 1, join = 'inner').sort_values(['CHR','POS'])
    out['P'] = sts.norm.sf(abs(out['Z'])) * 2
    out['BETA'] = out['Z'] / out['N'] / (out['AF1'] * (1-out['AF1'])) ** 0.5
    out['SE'] = out['BETA'] / out['Z']
    out.to_csv(args.out, sep = '\t', index = True, header = True)
    log.log(f'Output written to {args.out}')
    log.log('Analysis finished')

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description = 'A wrapper script to conduct genomic PCA or vanilla nw-GWAMA meta-analysis')
    parser.add_argument('pheno', nargs = '+', help = 'Phenotype groups to include in the analysis')
    parser.add_argument('-i','--in', dest = '_in', default = '../gwa', help = 'Directory of input GWAS summary statistics')
    parser.add_argument('--pca', action = 'store_true', help = 'Use PCA to estimate weights')
    parser.add_argument('--nw', action = 'store_true', help = 'Use n-weighted meta-analysis to estimate weights')
    parser.add_argument('--low_memory', action = 'store_true', help = 'Use low-memory mode')
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