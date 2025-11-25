#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2024-11-12
Version 2: 2025-08-22

Conducts MAGMA for single GWAS summary statistics

Requires following inputs: 
    GWAS summary statistics (single file)
    MAGMA annotation files (including H-MAGMA)
    MAGMA binary
'''

def main(args):
    import os
    from fnmatch import fnmatch
    
    # array submitter
    from _utils.slurm import array_submitter
    submitter = array_submitter(name = 'annot_magma_'+'_'.join(args.pheno),timeout = 360)
    if not os.path.isdir(args.out): os.mkdir(args.out)
    
    # find phenotypes
    from _utils.path import find_gwas
    pheno = find_gwas(args.pheno, long = True)

    # find annotations
    annots = [(f'{args.annot}/{x}', x.replace('.genes.annot','')) for x in os.listdir(args.annot) if x.endswith('.genes.annot')]

    for g, p in pheno:
        outdir = f'{args.out}/{g}'
        if not os.path.isdir(outdir): os.system(f'mkdir -p {outdir}')

        for annot, annot_prefix in annots:
            out_prefix = f'{outdir}/{p}.{annot_prefix}'
            if os.path.isfile(f'{out_prefix}.genes.out') and not args.force: continue
            submitter.add(f'{args.magma} --bfile {args.bfile} --gene-annot {annot} '+
                f'--pval {args._in}/{g}/{p}.fastGWA ncol=N --out {out_prefix}')
    submitter.submit()
    return submitter

if __name__ == '__main__':
    from _utils.slurm import slurm_parser
    parser = slurm_parser(
      description = 'This file batch runs the MAGMA pipeline for all IDPs')
    parser.add_argument('pheno', nargs = '*', help = 'Phenotypes')
    parser.add_argument('-i','--in', dest = '_in', help = 'Directory containing all phenotypes',
      default = '../gwa/')
    parser.add_argument('-a', '--annot', dest ='annot', help = 'directory to annotation files',
      default = '../toolbox/hmagma')
    parser.add_argument('-b', '--bfile', dest = 'bfile', help = 'bed binary to use in magma',
      default = '../toolbox/magma/g1000_eur')
    parser.add_argument('-o', '--out', dest = 'out', help = 'output directory',
      default = '../annot/magma')
    parser.add_argument('--magma', dest = 'magma', help = 'MAGMA executable',
      default = '../toolbox/magma/magma')
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
      default = False, action = 'store_true')
    args = parser.parse_args()
    
    # path normalisation
    import os
    for arg in ['_in','annot','bfile','out','magma']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    
    from _utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheng/%pheno_%maf.fastGWA', __file__)
    proj.add_output(args.out+'/%pheng/%pheno_%maf.log',__file__)
    try: main(args)
    except: cmdhistory.errlog()