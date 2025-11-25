#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-11-02

Conducts down-stream interval-based enrichment for LAVA results

Requires following inputs: 
    LAVA output
    inrich default resources
Outputs:
    Main analysis (pathway enrichments)
    Interval-gene-target analysis
'''

import os
import warnings

from matplotlib.pylab import f

def main(args):
    from ._utils.slurm import array_submitter
    submitter = array_submitter(name = 'gcorr_lava_inrich_'+'_'.join(args.p1)+'_'+'_'.join(args.p2),
        n_cpu = 1, timeout = 60, partition = 'sapphire')
    
    from ._utils.path import find_gwas, pair_gwas
    exposures = find_gwas(args.p1, dirname = args.gwa, ext = 'fastGWA', clump = True)
    outcomes = find_gwas(args.p2, dirname = args.gwa, ext = 'fastGWA', clump = True)
    pair = pair_gwas(exposures, outcomes)
    
    for g1, p1s, g2, p2s in pair:
        if g1 > g2: g1, p1s, g2, p2s = g2, p2s, g1, p1s
        for p1 in p1s:
            lava_file = f'{args._in}/{g1}.{g2}/{g1}_{p1}.{g2}.lava.txt'
            if not os.path.isfile(lava_file):
                warnings.warn(f'Missing LAVA output for {g1}/{p1} and {g2}, please run:\n' +
                    f'gcorr_lava_batch.py -p1 {g1} -p2 {g2} -i {os.path.dirname(args._in)}')
                continue
            out_prefix = f'{args.out}/{g1}.{g2}/{g1}_{p1}.{g2}.lava_inrich'
            os.makedirs(os.path.dirname(out_prefix), exist_ok = True)
            if (not args.force) and os.path.isfile(f'{out_prefix}.txt'): continue
            cmd = f'python gcorr_lava_inrich.py -i {lava_file} --inrich {args.inrich} -o {out_prefix}'
            submitter.add(cmd)
    submitter.submit()

if __name__ == '__main__':
    from ._utils.slurm import slurm_parser
    parser = slurm_parser(description = 'This scripts investigates interval-based enrichments for LAVA-identified genomic regions')
    parser.add_argument('-p1', nargs = '+', required = True, help = 'Exposures')
    parser.add_argument('-p2', nargs = '*', default = [], help = 'Outcomes, leave blank to run pairwise correlations across exposures')
    parser.add_argument('-i', '--in', dest = '_in', help = 'LAVA output directory',
        default = '../gcorr/lava/')
    parser.add_argument('-g','--gwa', help = 'Directory containing GWAS summary statistics',
        default = '../gwa/')
    parser.add_argument('--inrich', 
        help = 'folder of the inrich binary and resources, should contain resources/genes.txt and resources/snps.txt',
        default = '/rds/project/rds-Nl99R8pHODQ/toolbox/inrich')
    parser.add_argument('-o', '--out', help = 'output directory',
        default = '../gcorr/lava_inrich/')
    parser.add_argument('-f','--force', action = 'store_true', help = 'Force overwrite')
    args = parser.parse_args()
    for key in ['_in','gwa','inrich','out']:
        setattr(args, key, os.path.realpath(getattr(args, key)))
    
    from ._utils import cmdhistory, logger
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()