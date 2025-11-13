#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2025-11-04

A wrapper script to conduct GSA-mixer analysis

Requires following inputs:
    LDSC-munged GWAS summary statistics
    Gene sets
    Default GSA-mixer resources under <mixer directory>/resources
Outputs:
    Univariate GSA-mixer (if p2 is not specified)
    Bivariate GSA-mixer (if p2 is specified)
'''

import os
import tempfile

def split_sumstats(input, ref, output):
    import pandas as pd

    chroms = list(range(1,23))
    df = pd.read_table(input)
    for chrom in chroms:
        ref_df = pd.read_table(ref.replace('@', str(chrom)), header = None, usecols = [1])
        df_sub = df.loc[df['SNP'].isin(ref_df[1]),:]
        df_sub.dropna().to_csv(output.replace('@', str(chrom)), sep = '\t', index = False)
    return

def main(args):
    # tempdir
    mixer_py = f'{args.mixer}/precimed/mixer_dev.py'
    tmpdir = tempfile.mkdtemp()
    tmpdir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/temp/gsa_mixer'; os.makedirs(tmpdir, exist_ok = True)

    # find GWAS sumstats
    from _utils.path import find_gwas, pair_gwas
    exposures = find_gwas(args.p1, dirname = args._in, ext = 'sumstats', long = True)
    outcomes = find_gwas(args.p2, dirname = args._in, ext = 'sumstats', long = True)
    pheno = find_gwas(args.p1 + args.p2, dirname = args._in, ext = 'sumstats', long = True)
    pheno_pair = pair_gwas(exposures, outcomes)

    # find gene sets
    gene_sets = [x[:-4] for x in os.listdir(args.test) if x[-4:] == '.txt']

    # array submitter
    from _utils.slurm import array_submitter
    fit1_submitter = array_submitter(name = 'sc_gsamixer_fit1_'+'_'.join(args.p1)+'_'+'_'.join(args.p2),
        env = args.mixer, n_cpu = 8, timeout = 720)
    fit2_submitter = array_submitter(name = 'sc_gsamixer_fit2_'+'_'.join(args.p1)+'_'+'_'.join(args.p2),
        env = args.mixer, n_cpu = 8, timeout = 720, dependency=fit1_submitter)
    test1_submitter = array_submitter(name = 'sc_gsamixer_test1_'+'_'.join(args.p1)+'_'+'_'.join(args.p2),
        env = args.mixer, n_cpu = 8, timeout = 720, dependency=fit1_submitter)
    test2_submitter = array_submitter(name = 'sc_gsamixer_test2_'+'_'.join(args.p1)+'_'+'_'.join(args.p2),
        env = args.mixer, n_cpu = 8, timeout = 720, dependency=fit2_submitter)
    plsa_submitter = array_submitter(name = 'sc_gsamixer_plsa_'+'_'.join(args.p1)+'_'+'_'.join(args.p2),
        env = args.mixer, n_cpu = 16, timeout = 720)
    
    # mixer setup
    common_flags = ['--seed','20251104','--exclude-ranges','MHC','--threads','16',
        '--bim-file', f'{args.mixer}/resources/ldsc/1000G_EUR_Phase3_plink/chr@.bim',
        # '--loadlib-file', f'{args.mixer}/resources/ldsc/1000G_EUR_Phase3_plink/chr@.bin',
        '--ld-file', f'{args.mixer}/resources/ldsc/1000G_EUR_Phase3_plink/chr@.ld',
        '--annot-file', f'{args.mixer}/resources/ldsc/1000G_EUR_Phase3_plink/baseline_v2.2_chr@.annot.gz'
        # '--bim-file', f'{args.mixer}/resources/ukb_EUR_qc/chr@.bim',
        # '--loadlib-file', f'{args.mixer}/resources/ukb_EUR_qc/chr@.bin', # needs to be generated manually
        # '--annot-file', f'{args.mixer}/resources/ukb_EUR_qc/chr@.annot.gz', # needs to be generated manually
        ]
    
    # univariate steps
    for g, p in pheno:
        p_file = f'{args._in}/{g}/{p}.sumstats'
        temp_g = f'{tmpdir}/{g}_{p}.chr@.sumstats'
        if not os.path.isfile(temp_g.replace('@','22')):
            split_sumstats(p_file, f'{args.mixer}/resources/ldsc/1000G_EUR_Phase3_plink/chr@.bim', temp_g)
        if not os.path.islink(temp_g): os.symlink(p_file, temp_g)
        
        for gset in gene_sets:
            out_prefix = f'{args.out}/{g}/{g}_{p}.{gset}'
            os.makedirs(os.path.dirname(out_prefix), exist_ok = True)

            # plsa baseline model
            cmd0 = [
                'python', mixer_py, 'plsa', '--gsa-base', '--trait1-file', temp_g, 
                '--use-complete-tag-indices',
                '--go-file', f'{args.mixer}/resources/gsa-mixer-baseline-annot_10mar2023.csv',
                '--out', f'{out_prefix}.plsa.baseline'
            ] + common_flags
            cmd1 = [
                'python',mixer_py, 'plsa', '--gsa-full', '--trait1-file', temp_g, 
                '--use-complete-tag-indices',
                '--go-file', f'{args.mixer}/resources/gsa-mixer-gene-annot_10mar2023.csv',
                '--go-file-test', f'{args.test}/{gset}.txt', # f'{args.mixer}/resources/gsa-mixer-hybridLOO-annot_10mar2023.csv'
                '--load-params', f'{out_prefix}.plsa.baseline.json',
                '--out', f'{out_prefix}.plsa.full'
            ] + common_flags

            cmd = []
            if args.force or not os.path.isfile(f'{out_prefix}.plsa.baseline.json'):
                cmd.append(' '.join(cmd0))
            if args.force or not os.path.isfile(f'{out_prefix}.plsa.full.json'):
                cmd.append(' '.join(cmd1))
            if cmd: plsa_submitter.add(*cmd) # ensure that the two commands are run in order

            # univariate fit
            cmd = [
                'python', mixer_py, 'fit1', '--trait1-file', temp_g,
                '--out', f'{out_prefix}.fit', '--go-file', f'{args.test}/{gset}.txt'
            ] + common_flags
            if args.force or not os.path.isfile(f'{out_prefix}.fit.json'):
                fit1_submitter.add(' '.join(cmd))
            
            # univariate test
            cmd = [
                'python', mixer_py, 'test1', '--trait1-file', temp_g,
                '--load-params', f'{out_prefix}.fit.json',
                '--out', f'{out_prefix}.test', '--go-file', f'{args.test}/{gset}.txt'
            ] + common_flags
            if args.force or not os.path.isfile(f'{out_prefix}.test.json'):
                test1_submitter.add(' '.join(cmd))

    for g1, p1, g2, p2 in pheno_pair:
        if (g1 > g2) or (g1 == g2 and p1 > p2):
            g1, g2, p1, p2 = g2, g1, p2, p1
            
        for gset in gene_sets:
            if g1 == g2 and p1 == p2: continue
            out_prefix = f'{args.out}/{g1}.{g2}/{g1}_{p1}.{g2}_{p2}.{gset}'
            os.makedirs(os.path.dirname(out_prefix), exist_ok = True)
            temp_g1 = f'{tmpdir}/{g1}_{p1}.chr@.sumstats'
            temp_g2 = f'{tmpdir}/{g2}_{p2}.chr@.sumstats'
            if not os.path.isfile(temp_g1.replace('@','22')): 
                raise FileNotFoundError(f'File not found: {temp_g1.replace("@","22")}, please check sumstats')
            if not os.path.isfile(temp_g2.replace('@','22')): 
                raise FileNotFoundError(f'File not found: {temp_g2.replace("@","22")}, please check sumstats')
            
            # bivariate fit
            cmd = [
                'python', mixer_py, 'fit2', '--trait1-file', temp_g1, '--trait2-file', temp_g2,
                '--trait1-params', f'{args.out}/{g1}/{g1}_{p1}.{gset}.fit.json',
                '--trait2-params', f'{args.out}/{g2}/{g2}_{p2}.{gset}.fit.json',
                '--out', f'{out_prefix}.fit', '--go-file', f'{args.test}/{gset}.txt'
            ] + common_flags
            if args.force or not os.path.isfile(f'{out_prefix}.fit.json'):
                fit2_submitter.add(' '.join(cmd))
            
            # bivariate test
            cmd = [
                'python', mixer_py, 'test2', '--trait1-file', temp_g1, '--trait2-file', temp_g2,
                '--load-params', f'{out_prefix}.fit.json',
                '--out', f'{out_prefix}.test', '--go-file', f'{args.test}/{gset}.txt'
            ] + common_flags
            if args.force or not os.path.isfile(f'{out_prefix}.test.json'):
                test2_submitter.add(' '.join(cmd))
    plsa_submitter.submit()
    fit1_submitter.submit()
    fit2_submitter.submit()
    test1_submitter.submit()
    test2_submitter.submit()
    
if __name__ == '__main__':
    from _utils.slurm import slurm_parser
    parser = slurm_parser(description = 'This script conducts GSA-mixer analysis between given pairs of phenotypes')
    parser.add_argument('-p1','--p1', nargs = '+', required=True, help='Phenotype group 1')
    parser.add_argument('-p2','--p2', nargs = '*', default = [], help='Phenotype group 2 (if not specified, p2 = p1)')
    parser.add_argument('-i','--in', dest = '_in', help = 'Input directory containing LDSC-munged GWAS summary statistics',
        default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/gcorr/ldsc_sumstats/')
    parser.add_argument('--test', help='Path to gene set file for testing',
        default = '/rds/project/rds-Nl99R8pHODQ/multiomics/gene_set/gene_set_4mixer')
    parser.add_argument('--mixer', help='Path to GSA-mixer installation directory',
        default = '/rds/project/rds-Nl99R8pHODQ/toolbox/mixer')
    parser.add_argument('-o','--out', help = 'Output directory',
        default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/sc/gsa_mixer/')
    parser.add_argument('-f','--force', action = 'store_true', help = 'Force overwrite')
    args = parser.parse_args()

    import os
    for arg in ['_in','out','test']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))

    from _utils import cmdhistory, logger
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()