#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2023-07-22

This is a script for single-trait fine-mapping using polyfun-susie

Required input:
    GWAS summary statistics (scans directory for files)
'''

def main(args):
    # make output directories
    import os
    if not os.path.isdir(args.out): os.mkdir(args.out)
    
    # array submitter
    from ._utils.slurm import array_submitter
    submitter = array_submitter(name = f'finemap_{args.pheno[0]}', n_cpu = 2,
        env = '/rds/project/rds-Nl99R8pHODQ/toolbox/polyfun', timeout = 120)
    
    from ._utils.path import find_gwas, find_clump
    pheno = find_gwas(args.pheno, dirname = args._in, clump = True, long = True)
    for g, p in pheno:
        try: clump,_ = find_clump(g, p, dirname = args.clump, pval = args.p)
        except: Warning(f'No clump file found for {g}/{p} - skipping'); continue

        outdir = f'{args.out}/{g}'
        if not os.path.isdir(outdir): os.system(f'mkdir -p {outdir}')

        out_file = f'{args.out}/{g}/{p}.finemap.summary'
        if os.path.isfile(out_file) and (not args.force): continue

        cmd = ['python', 'finemap_by_trait.py', '-i', f'{args._in}/{g}/{p}.fastGWA', '-c', clump, 
            '-o', out_file, '-b', args.bfile, '-p', f'{args.p:.4e}', '--polyfun', args.polyfun]
        if args.force: cmd.append('-f')
        submitter.add(' '.join(cmd))
    submitter.submit()
    return submitter
    
if __name__ == '__main__':
    from ._utils.slurm import slurm_parser
    parser = slurm_parser(
      description = 'This programme batch runs the fine-map pipeline')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
    parser.add_argument('-i','--in', dest = '_in', help = 'Directory containing all summary stats',
      default = '../gwa/')
    parser.add_argument('-c','--clump', dest = 'clump', help = 'Directory containing all clump outputs',
      default = '../clump/')
    parser.add_argument('-o', '--out', dest = 'out', help = 'output directory',
      default = '../finemap/')
    parser.add_argument('-b', '--bfile', dest = 'bfile', help = 'bed binary',
      default = '../params/bed/')
    parser.add_argument('--polyfun', help = 'directory of POLYFUN tool',
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/polyfun/')
    parser.add_argument('-p', dest = 'p', help = 'p-value', default = 5e-8, type = float)
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
      default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','clump','bfile']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    
    from ._utils import path, cmdhistory, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheng/%pheno_%maf.fastGWA', __file__)
    proj.add_output(args.out+'/%pheng/%pheno_%maf.finemap.summary',__file__)
    try: main(args)
    except: cmdhistory.errlog()