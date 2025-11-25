#!/usr/bin/env python3
'''
Computes correlational models between IDPs and PGS
'''

def main(args):
    import os
    from fnmatch import fnmatch
    # from sklearn.linear_model import LinearRegression
    
    # array submitter
    from _utils.slurm import array_submitter
    submitter = array_submitter(
        name = f'prs_corr_{args.pheno[0]}', n_cpu = 1,
        timeout = 30,
        # debug = True
        )
    
    # input parsing
    f = ' -f' if args.force else ''
    if not os.path.isdir(args.out): os.mkdir(args.out)
    diagdir = args.out + '/diagnostics/'
    if not os.path.isdir(diagdir): os.mkdir(diagdir)
    
    flist = []
    for p in args.pheno:
        if os.path.isfile(f'{args._in}/{p}'):
            flist.append(p)
            continue
        if os.path.isfile(f'{args._in}/{p}.txt'):
            flist.append(f'{p}.txt')
            continue
        
        for x in os.listdir(args._in):
            if fnmatch(x, f'*{p}*.txt'):
                flist.append(x)
    print('Files to process:')
    for x in flist: print(x)
    
    # Correlation
    for x in flist:
      os.system('python '+
        f'prs_corr.py -i {args._in}/{x} -p {" ".join(args.prs)} --prsdir {args.prsdir} --dcov {args.dcov} --qcov {args.qcov} -o {args.out} {f}')
    submitter.submit()

if __name__ == '__main__':
    from _utils.slurm import slurm_parser
    parser = slurm_parser(description='Computes correlational plots between IDPs and PGS (batch)')
    parser.add_argument('pheno', nargs = '*', help = 'Phenotype files to process')
    parser.add_argument('-i','--in', dest = '_in', help = 'Input directory',
      default = '../pheno/ukb/')
    parser.add_argument('-p','--prs', dest = 'prs', help = 'PRS to select', nargs = '*',
      default = ['disorders','disorders_subtypes'])
    parser.add_argument('--prsdir', dest = 'prsdir', help = 'PRS score directory',
      default = '../prs/prs_score/')
    parser.add_argument('--dcov',dest = 'dcov', help = 'DISCRETE covariance file',
      default = '../params/ukb_dcov.txt')
    parser.add_argument('--qcov',dest = 'qcov', help = 'QUANTITATIVE covariance file',
      default = '../params/ukb_qcov.txt')
    parser.add_argument('-o','--out', dest = 'out', help = 'Output directory',
      default = '../prs/prs_corr/')
    parser.add_argument('-f','--force', dest = 'force', help = 'force overwrite',
                        action = 'store_true', default = False)
    args = parser.parse_args()
    import os
    for arg in ['_in','out','prsdir','dcov','qcov']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    args.pheno.sort()
    args.prs.sort()
    
    from _utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheng.txt', __file__)
    proj.add_output(args.out+'/%pheno_.*.txt', __file__)
    try: main(args)
    except: cmdhistory.errlog()