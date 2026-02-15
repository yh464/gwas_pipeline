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

from hashlib import sha256
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

def sep_chr(input_args):
    g, p, in_dir, tmpdir = input_args
    if all([os.path.exists(f'{tmpdir}/{chrom}/{g}/{p}.parquet') for chrom in range(1,23)]):
        log.log(f'{g}/{p} already separated into chromosomes')
        if os.path.exists(f'{tmpdir}/23/{g}/{p}.parquet'): return list(range(1,24))
        elif os.path.exists(f'{tmpdir}/X/{g}/{p}.parquet'): return list(range(1,23)) + ['X']
        else: return list(range(1,23))
        
    df = pd.read_table(f'{in_dir}/{g}/{p}.fastGWA', usecols = 
        lambda x: x.upper() in (['CHR','SNP','POS','A1','A2','BETA','OR','SE','N','AF1']),
        index_col = ['SNP'], dtype = {
            'CHR': 'category', 'POS': np.int32, 'SNP': str, 'A1': 'category', 'A2': 'category',
            'BETA': np.float32, 'OR': np.float32, 'SE': np.float32, 'N': np.float32, 'AF1': np.float32
    })
    for chrom, df_chr in df.groupby('CHR'):
        os.makedirs(f'{tmpdir}/{chrom}/{g}', exist_ok = True)
        df_chr.to_parquet(f'{tmpdir}/{chrom}/{g}/{p}.parquet', index = True)
    gc.collect()
    return df['CHR'].unique().tolist()

def read_sumstats(input_args):
    g, p, in_dir, weight, = input_args
    if os.path.isfile(f'{in_dir}/{g}/{p}.fastGWA'):
        df = pd.read_table(f'{in_dir}/{g}/{p}.fastGWA', usecols = 
            lambda x: x.upper() in (['CHR','SNP','POS','A1','A2','BETA','OR','SE','N','AF1']),
            index_col = ['SNP'], dtype = {
                'CHR': 'category', 'POS': np.int32, 'SNP': str, 'A1': 'category', 'A2': 'category',
                'BETA': np.float32, 'OR': np.float32, 'SE': np.float32, 'N': np.float32, 'AF1': np.float32
        })
    elif os.path.isfile(f'{in_dir}/{g}/{p}.parquet'):
        df = pd.read_parquet(f'{in_dir}/{g}/{p}.parquet')
        df.index.name = 'SNP'
    else: 
        log.warn(f'File not found for {in_dir}/{g}/{p}, returning empty dataframes')
        return (pd.Series(dtype = np.float32), 
                  pd.Series(dtype = np.float32), 
                  pd.Series(dtype = np.float32), 
                  pd.Series(dtype = np.float32), 
                  pd.DataFrame(columns = ['CHR','POS','A1','A2'], index = [], dtype = {
                      'CHR': 'category', 'POS': np.int32, 'A1': 'category', 'A2': 'category'}))
    df = df.loc[~df.index.duplicated(keep = False), :].sort_index()
    df.columns = [x.upper() for x in df.columns]
    if 'OR' in df.columns and not 'BETA' in df.columns: df['BETA'] = np.log(df['OR'])
    n = df['N']
    w = (weight * (n ** 0.5))
    z = (df['BETA'] * w / df['SE'])
    af = (df['AF1'] * n)
    snpinfo = df[['CHR','POS','A1','A2']].copy()
    del df
    gc.collect()
    return w, z, af, n, snpinfo

def gwama(pheno, weight_list, z_list, af_list, n_list, snpinfo_list, gcovint):
    log.log('Merging variant information across all summary statistics')
    snpinfo = pd.concat(snpinfo_list, axis = 0).drop_duplicates().sort_index()
    if snpinfo.duplicated(['CHR','POS']).any():
        log.warn(f'{snpinfo.duplicated(["CHR","POS"]).sum()} variants have different alleles across files')
        out_missnp = snpinfo.loc[snpinfo.duplicated(['CHR','POS'], keep = False), :]
        snpinfo = snpinfo.drop_duplicates(['CHR','POS'], keep = False)
        log.log(f'{snpinfo.shape[0]} variants with consistent alleles across files retained for meta-analysis')
    else: out_missnp = pd.DataFrame(columns = snpinfo.columns)
    del snpinfo_list
    gc.collect()

    log.log('Calculating meta-analytic Z-scores')
    n_total = pd.concat(list(n_list), axis = 1, ignore_index = True).fillna(0).sum(axis = 1).rename('N')
    out_z = pd.concat(list(z_list), axis = 1, ignore_index = True).fillna(0).sum(axis = 1).rename('Z')
    af1 = (pd.concat(list(af_list), axis = 1, ignore_index = True).fillna(0).sum(axis = 1) / n_total).rename('AF1')
    weight = pd.concat(list(weight_list), axis = 1, ignore_index = True).fillna(0)
    weight.columns = pd.MultiIndex.from_tuples(pheno, names = ['group','pheno'])
    gcovint = gcovint.loc[weight.columns, weight.columns].fillna(0).values
    div_coef = np.einsum('ij, jk, ik -> i', weight.values, gcovint, weight.values, optimize = 'optimal') ** 0.5
    out_z = out_z / div_coef
    out = pd.concat([snpinfo, af1, out_z, n_total], axis = 1, join = 'inner').sort_values(['CHR','POS'])
    out['P'] = sts.norm.sf(abs(out['Z'])) * 2
    out['BETA'] = out['Z'] / out['N'] / (out['AF1'] * (1-out['AF1'])) ** 0.5
    out['SE'] = out['BETA'] / out['Z']
    return out, out_missnp

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
    pool = Pool(min(args.threads, len(pheno)))
    log.log(f'Starting parallel pool using {min(args.threads, len(pheno))} threads')

    # read SNP info for each trait
    if args.sep_chr:
        tmpdir = os.path.realpath('../temp/gwama_sep_chr')
        parallel_args = [(g, p, args._in, tmpdir) for g, p in pheno]
        chroms = list(tqdm(pool.imap(sep_chr, parallel_args, chunksize = min(args.threads, len(parallel_args))), 
            total = len(parallel_args), 
            desc = 'Separating chromosomes for each trait'))
        chroms = list(set([chrom for sublist in chroms for chrom in sublist]))

        out = []; out_missnp = []
        for chrom in tqdm(chroms, desc = 'Processing each chromosome'):
            parallel_args = [(g, p, f'{tmpdir}/{chrom}', weights.loc[(g,p)]) for g, p in pheno]
            chr_temp = f'{tmpdir}/{chrom}/{sha256(str(pheno).encode()).hexdigest()[:12]}'
            if os.path.isfile(f'{chr_temp}.ss.parquet') and os.path.isfile(f'{chr_temp}.missnp.parquet'):
                log.log(f'Chromosome {chrom}: found existing processed files, loading from disk')
                out_chr = pd.read_parquet(f'{chr_temp}.ss.parquet')
                missnp_chr = pd.read_parquet(f'{chr_temp}.missnp.parquet')
                out.append(out_chr)
                out_missnp.append(missnp_chr)
            else:
                chr_out = list(tqdm(pool.imap(read_sumstats, parallel_args, chunksize = min(args.threads, len(parallel_args))), 
                    total = len(parallel_args), 
                    desc = f'Processing chromosome {chrom}'))
                weight_list, z_list, af_list, n_list, snpinfo_list = zip(*chr_out)
                out_chr, missnp_chr = gwama(pheno, weight_list, z_list, af_list, n_list, snpinfo_list, gcovint)
                
                out.append(out_chr)
                out_missnp.append(missnp_chr)
                out_chr.to_parquet(f'{chr_temp}.ss.parquet')
                missnp_chr.to_parquet(f'{chr_temp}.missnp.parquet')
                del weight_list, z_list, af_list, n_list, snpinfo_list, chr_out
                gc.collect()
        out = pd.concat(out, axis = 0).sort_values(['CHR','POS'])
        out_missnp = pd.concat(out_missnp, axis = 0)

    else:
        parallel_args = [(g, p, args._in, weights.loc[(g,p)]) for g, p in pheno]
        out = list(tqdm(pool.imap(read_sumstats, parallel_args, chunksize = min(args.threads, len(parallel_args))), 
            total = len(parallel_args), 
            desc = 'Reading summary statistics and aggregating weighted Z-scores'))
        weight_list, z_list, af_list, n_list, snpinfo_list = zip(*out)
        out, out_missnp = gwama(pheno, weight_list, z_list, af_list, n_list, snpinfo_list, gcovint)
    out.to_csv(args.out, sep = '\t', index = True, header = True)
    if out_missnp.shape[0] > 0:
        log.warn(f'{out_missnp.shape[0]} variants have different alleles across files, written to {args.out.replace(".fastGWA", ".missnp")}')
        out_missnp.to_csv(args.out.replace('.fastGWA', '.missnp'), sep = '\t', index = False)
    log.log(f'Output written to {args.out}')
    log.log('Analysis finished')

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description = 'A wrapper script to conduct genomic PCA or vanilla nw-GWAMA meta-analysis')
    parser.add_argument('pheno', nargs = '+', help = 'Phenotype groups to include in the analysis')
    parser.add_argument('-i','--in', dest = '_in', default = '../gwa', help = 'Directory of input GWAS summary statistics')
    parser.add_argument('--pca', action = 'store_true', help = 'Use PCA to estimate weights')
    parser.add_argument('--nw', action = 'store_true', help = 'Use n-weighted meta-analysis to estimate weights')
    parser.add_argument('--sep_chr', action = 'store_true', help = 'Separate chromosomes, useful for low memory mode')
    parser.add_argument('--threads', type = int, default = 16, help = 'Number of threads to use for parallel processing')
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