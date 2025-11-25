#!/usr/bin/env python3
'''
generates manhattan plots for all fastGWA files in a directory
'''

def main(args):
    import os
    from fnmatch import fnmatch
    from _utils.slurm import array_submitter
    
    f = f'-p {args.pval}'
    if args.force: f += ' -f'
    if args.a: f += ' -a'
    if len(args.pheno) == 0:
        args.pheno = os.listdir(args._in)
    if not os.path.isdir(args.out): os.system(f'mkdir -p {args.out}')
    
    submitter = array_submitter(
        name = f'gwa_manhattan_{args.pheno[0]}',
        partition = 'icelake-himem', # must use himem, using icelake will raise out-of-memory error [125]
        timeout = 60,
        debug = False
        )
    
    for y in args.pheno:
        os.chdir(args._in)
        if not os.path.isdir(f'{args.out}/{y}'): os.mkdir(f'{args.out}/{y}')
        if not os.path.isdir(y): continue # sanity check especially for default y = ls(args._in)
        os.chdir(y)
        flist = []
        for x in os.listdir():
            if fnmatch(x, '*.fastGWA') and not fnmatch(x, '*X.fastGWA'): 
                flist.append(x)
        
        if not os.path.isdir(args.out+y):
            os.system(f'mkdir -p {args.out}/{y}')
        
        for x in flist:
            out_fname = f'{args.out}/{x}'.replace('.fastGWA','.manhattan.png')
            if os.path.isfile(out_fname) and not args.force: continue
            submitter.add('python '+
              f'gwa_manhattan.py {y} --file {x} -i {args._in} -o {args.out} {f}')
    
    submitter.submit()
    # submitter.debug()
        
if __name__ == '__main__':
    from _utils.slurm import slurm_parser
    parser = slurm_parser(description = 
      'This programme compiles Manhattan plots for all fastGWA output files for a single phenotype file')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
    parser.add_argument('-i','--in', dest = '_in', help = 'GWA file directory',
      default = '../gwa/')
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
      default = '../gwa/manhattan/')
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
      default = False, action = 'store_true')
    parser.add_argument('-p','--pval', help = 'p-value threshold', type = float, default = 5e-8)
    parser.add_argument('-a','--autosome-only',dest = 'a', help = 'exclude sex chromosomes',
      default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['_in','out']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in, __file__)
    proj.add_output(args.out, __file__)
    try: main(args)
    except: cmdhistory.errlog()