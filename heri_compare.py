#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-02-25

Compares LDSC and GREML based heritability
'''

def main(args):
    import pandas as pd
    from fnmatch import fnmatch
    import numpy as np
    import scipy.stats as sts
    from logparser import parse_h2_log, parse_greml_h2_log
    
    summary = []
    for p in args.pheno:
        for x in os.listdir(f'{args.ldsc}/{p}'):
            if not fnmatch(x,'*.sumstats'): continue
            prefix = x.replace('.sumstats','').replace('.gz','')
            ldsc = f'{args.ldsc}/{p}/{prefix}.h2.log'
            ldsc_h2, ldsc_se = parse_h2_log(ldsc)
            greml = f'{args.greml}/{p}/{prefix}.greml.hsq'
            greml_h2, greml_se = parse_greml_h2_log(greml)
            
            summary.append(pd.DataFrame(dict(
                group = [p], phenotype = prefix,
                ldsc_h2 = ldsc_h2, ldsc_se = ldsc_se, ldsc_p = 1-sts.chi2.cdf(ldsc_h2**2/ldsc_se**2, df = 1),
                greml_h2 = greml_h2, greml_se = greml_se, greml_p = 1-sts.chi2.cdf(greml_h2**2/greml_se**2, df = 1),
                )))
    summary = pd.concat(summary)
    summary.insert(loc = 5, column = 'ldsc_q', value = sts.false_discovery_control(summary.ldsc_p))
    summary['greml_q'] = sts.false_discovery_control(summary.greml_p)
    from _utils.path import normaliser
    norm = normaliser()
    norm.normalise(summary).to_csv(f'{args.out}/h2_se_'+'_'.join(args.pheno) + '.txt', sep = '\t', index = False)
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'This programme generates a summary table for LDSC and GREML heritability estimates')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
    parser.add_argument('--greml', dest = 'greml', help = 'GREML logfile directory',
        default = '../gwa/')
    parser.add_argument('--ldsc', help = 'LDSC logfile directory',
        default = '../gcorr/ldsc_sumstats/')
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
        default = '../gcorr/')
    args = parser.parse_args()
    args.pheno.sort()
    import os
    for arg in ['greml','ldsc','out']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    
    from _utils import cmdhistory, logger
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()