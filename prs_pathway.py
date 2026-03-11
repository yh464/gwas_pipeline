#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2026-03-09

Python wrapper to run pathway-specific PRS by partitioning PRS-cs effect sizes

Inputs:
    PRS effect sizes from PRS-cs (columns: CHR, SNP, POS, A1, A2, BETA), separated by chromosome
    PLINK .bed files of target population
    A gene set / gene score file
'''

import os
import pandas as pd
import numpy as np
from tqdm import tqdm
from _utils.logger import logger
log = logger()

def main(args):
    from bed_reader import open_bed # pyright: ignore[reportMissingImports]
    from _utils.genetools import ensg_to_loc

    geno = open_bed(f'{args.bed}.bed', properties = {'father': None, 'mother': None, 'sex': None, 'pheno': None})
    bim = pd.DataFrame(dict(
        CHR = geno.chromosome, POS = geno.bp_position, A1 = geno.allele_1, A2 = geno.allele_2, SNP = geno.sid
    )).set_index('SNP')
    fam = pd.DataFrame(dict(
        FID = geno.fid, IID = geno.iid
    )).set_index('IID')

    # read the gene set / gene score file
    hdr = open(args._in).readline().strip().split('\t')[1]
    if hdr.startswith('ENSG'): # gene_set gene1 gene2 ... format
        gsets = pd.concat([
            pd.DataFrame(columns = [x.split()[0]], index = x.split()[1:], data = 1) for x in open(args._in).read().splitlines() if x.split()[0] == args.gset
        ], axis = 1).fillna(0)
    else: gsets = pd.read_table(args._in, index_col = 0, usecols = [0, args.gset]).fillna(0) # index_col = gene, columns = gene sets, values = gene scores
    if gsets.shape[1] != 1: log.error(f'Multiple gene sets found in the input file. Please specify one gene set using --gset. Found gene sets: {gsets.columns.tolist()}')
    gsets = pd.concat(gsets).replace(0, pd.NA).dropna() # remove all-zero SNPs and fill NAs with zero
    log.log(f'{gsets.shape[0]} genes were found in the gene set {args.gset}')

    gene_loc = ensg_to_loc(gsets.index, build = args.build, window_up = args.window[0], window_down = args.window[1])

    # read PRS-cs effect sizes
    log.log('Reading PRS-cs effect sizes')
    effsize = pd.concat([
        pd.read_table(args.prs.replace('chr*', f'chr{chrom}'), header = None, names = ['CHR','SNP','POS','A1','A2','BETA']).set_index('SNP')
        # use lower case to be compatible with pyplink API
        for chrom in range(1, 23)
    ])

    # harmonise allele orders for effect sizes and target genotypes
    log.log('Harmonising allele orders for effect sizes and target genotypes')
    effsize = effsize.loc[effsize.index.intersection(bim.index),:]
    effsize = pd.concat([effsize, bim.loc[effsize.index,:][['A1','A2']].rename(columns = {'A1': 'A1_bim', 'A2': 'A2_bim'})], axis = 1)
    flip = (effsize['A1'] == effsize['A2_bim']) & (effsize['A2'] == effsize['A1_bim'])
    keep = (effsize['A1'] == effsize['A1_bim']) & (effsize['A2'] == effsize['A2_bim'])
    effsize = effsize[flip | keep].copy()
    effsize.loc[flip, 'BETA'] *= -1
    effsize = effsize.drop(columns = ['A1','A2']).rename(columns = {'A1_bim': 'A1', 'A2_bim': 'A2'})

    # find all SNPs within each gene to get SNP-level weights
    log.log('Generating SNP-level weights for each gene x pathway')
    snp_weights = []
    for gene, scores in gsets.iterrows():
        chrom = gene_loc.loc[gene,'CHR']; start = gene_loc.loc[gene,'START']; stop = gene_loc.loc[gene,'STOP']
        snps = bim.index[(bim['CHR'] == chrom) & (bim['POS'] >= start) & (bim['POS'] <= stop)]
        snp_weights.append(pd.DataFrame(index = snps, columns = gsets.columns, data = scores.values))
    n_snp = snp_weights.shape[0]
    log.log(f'{n_snp} SNPs were mapped to the gene set')

    # generate 10000 permuted SNP sets for permuted pathway PRS
    rng = np.random.default_rng(seed = 19260817)
    permuted_sets = [rng.choice(bim['SNP'], size = n_snp, replace = False) for _ in range(10000)]

    # generate pathway-specific PRS using snp-level weights and the null distribution using the permuted SNP sets
    log.log('Generating pathway-specific PRS and null distribution')
    out = pd.DataFrame(index = fam.index, columns = ['score_norm'] + [f'score_null_{i}' for i in range(10000)])

    for col, snpset in tqdm(zip(out.columns, [snp_weights.index] + permuted_sets), desc = 'Generating PRS'):
        snpidx = bim.index.get_indexer(snpset)
        geno_set = geno.read(snpidx) # columns are SNPs, rows are individuals
        geno_set = geno_set.nan_to_num(np.nanmean(geno_set, axis = 0)) # AF1 imputation for missing genotypes
        tmp = geno_set @ (snp_weights.values * effsize.loc[snpset,'BETA'].values) # SNP-level weight = PRS-cs effect size * gene set/score
        del geno_set
        out[col] = tmp / tmp.std() # standardise the PRS to variance = 1
    
    out.insert(0, 'FID', fam['FID'])
    out.insert(1, 'IID', fam.index)
    out.to_csv(args.out, index = False, sep = '\t')