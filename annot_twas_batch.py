#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2026-05-06

Conducts TWAS using FUSION for GWAS summary statistics

Requires following inputs: 
    GWAS summary statistics (LDSC formatted)
    FUSION reference genome (LDREF)
    FUSION weights directory and names
'''

from _utils import cmdhistory, path, logger
import os
log = logger.logger()
proj = path.project()

def main(args):
    from _utils.path import find_gwas
    from _utils.slurm import array_submitter
    import pandas as pd
    pheno = find_gwas(args.pheno, long = True)
    submitter = array_submitter(name = 'annot_twas_'+'_'.join([x[1] for x in pheno]), 
        env = 'gentoolsr', wd = args.fusion, timeout = 120, n_cpu = 4, partition = 'sapphire')
    proj.register('annot_twas', 'annot/twas/$group/$pheno/$weight.twas.txt')
    
    # select weights to be used
    l = [f.replace('.pos','') for f in os.listdir(args.weights) if f.endswith('.pos')]
    log.log('Following TWAS weights are found in the weights directory:')
    for i, f in enumerate(l):
        log.log(f'    {i+1}: {f}')
    default_weights = [x for x in ['cortex.nofilter', 'cross_disorder', 'sCCA1'] if x in l]
    selected_weights = [l[int(x)-1].replace('.pos', '') for x in input('Please select weights to be used, separated by whitespace:\n' + \
        str(default_weights) + '\n')] or default_weights
    if not selected_weights: log.error('No weights selected, exiting')
    print()
    log.log('Selected following weights:')
    for w in selected_weights: log.log(f'    {w}')

    for g, p in pheno:
      sumstats_file = proj.to_pathname('ldsc_sumstats', g, p)
      for weight in selected_weights:
        out_filename = proj.to_pathname('annot_twas', g, p, weight = weight)
        outdir = os.path.dirname(out_filename)
        os.makedirs(f'{outdir}/{weight}', exist_ok = True)

        chrom_files = []
        for chrom in range(1, 23): # autosomes only
          chrom_file = f'{outdir}/{weight}/chr{chrom}.twas.txt'
          if os.path.isfile(chrom_file) and not args.force: chrom_files.append(chrom_file); continue
          cmd = ['Rscript', f'{args.fusion}/FUSION.assoc_test.R', '--sumstats', sumstats_file, 
              '--weights', f'{args.weights}/{weight}.pos', '--weights_dir', args.weights, 
              '--ref_ld_chr', args.ldref, '--chr', str(chrom), '--out', chrom_file]
          submitter.add(' '.join(cmd))
        
        if len(chrom_files) == 22: pd.concat([pd.read_table(f) for f in chrom_files], axis = 0).to_csv(out_filename, sep = '\t', index = False)
    submitter.submit()
    return submitter

if __name__ == '__main__':
    from _utils.slurm import slurm_parser
    parser = slurm_parser(
        description = 'This file batch runs the TWAS pipeline for all IDPs')
    parser.add_argument('pheno', nargs = '*', help = 'Phenotypes')
    # input and output directories are fixed because these are important intermediate files
    parser.add_argument('--fusion', help = 'Directory to FUSION scripts', default = '../toolbox/fusion')
    parser.add_argument('--ldref', help = 'FUSION LD reference, PLINK format') # default to None, as relative to the fusion directory
    parser.add_argument('--weights', help = 'FUSION weights directory') # default to None, as relative to the fusion directory
    parser.add_argument('--force', dest = 'force', help = 'force output',
        default = False, action = 'store_true')
    args = parser.parse_args()

    if args.ldref is None: args.ldref = os.path.join(args.fusion, 'ldref', 'chr')
    if args.weights is None: args.weights = os.path.join(args.fusion, 'weights')
    for arg in ['fusion', 'ldref', 'weights']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))

    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()