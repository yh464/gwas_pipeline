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
    from ._utils.slurm import array_submitter
    submitter = array_submitter(
        name = 'local_spin', n_cpu = 1,
        timeout = 30, env = 'rwd',
        debug = True
        )
    
    force = '-f' if args.force else ''
    
    for x in args._in:
        x = os.path.realpath(x)
        d = os.path.dirname(x)
        p = os.path.basename(x)
        for f in os.listdir(d):
            if fnmatch(f,p) and not fnmatch(f,'*yeo.*') and not fnmatch(f,'*mes.*'):
                lr_out = f'{d}/{f}'.replace('.txt','').replace('.csv','')+'.spinlr.yeo.txt'
                out = lr_out.replace('spinlr','spin')
                if not os.path.isfile(lr_out) or args.force:
                    submitter.add(f'Rscript local_spin_lr.r -i {d}/{f} {force}')
                if not os.path.isfile(out) or args.force:
                    submitter.add(f'Rscript local_spin.r -i {d}/{f} {force}')
    
    submitter.submit()
    
if __name__ == '__main__':
    from ._utils.slurm import slurm_parser
    parser = slurm_parser(
        description = 'This programme batch runs spin permutations')
    parser.add_argument('-i','--in', dest = '_in', help = 'input files, can have wildcards, nargs=*',
        nargs = '*', 
        # default = ['../local_corr/*20??_rg.csv',
        #            '../local_corr/local_h2_summary.csv',
        #            '../local_corr/global_rg.csv',
        #            '../local_corr/asym_rg.csv',
        #            '../prs/prs_corr/*20??_beta.csv',
        #            '../clump/*3e-11_overlaps.txt']
        default = ['../brainplots/*.txt','../brainplots/*.csv']
        )
    parser.add_argument('-f','--force',dest = 'force', help = 'force overwrite',
        default = False, action = 'store_true')
    args = parser.parse_args()
    
    from ._utils import cmdhistory
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()