#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2024-10-29

This scripts batch submits files for spin permutation test using the R script
'''

def main(args):
    import os
    from fnmatch import fnmatch
    
    # array submitter
    from _utils import array_submitter
    submitter = array_submitter.array_submitter(
        name = 'local_spin', n_cpu = 1,
        timeout = 10, env = 'rwd',
        debug = False
        )
    
    force = '-f' if args.force else ''
    
    for x in args._in:
        x = os.path.realpath(x)
        d = os.path.dirname(x)
        p = os.path.basename(x)
        for f in os.listdir(d):
            if fnmatch(f,p):
                submitter.add(f'Rscript local_spin_lr.r -i {d}/{f} {force}')
                submitter.add(f'Rscript local_spin.r -i {d}/{f} {force}')
    
    submitter.submit()
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description = 'This programme batch runs spin permutations')
    parser.add_argument('-i','--in', dest = '_in', help = 'input files, can have wildcards, nargs=*',
        nargs = '*', 
        default = ['../local_corr/*20??_rg.csv',
                   '../prs/prs_corr/*20??_beta.csv',
                   '../clump/*3e-11_overlaps.txt'])
    parser.add_argument('-f','--force',dest = 'force', help = 'force overwrite',
        default = False, action = 'store_true')
    args = parser.parse_args()
    
    from _utils import cmdhistory
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()