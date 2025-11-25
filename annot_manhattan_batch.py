#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-02-26
Summarises gene-set level enrichment for HMAGMA and MAGMA outputs

Preceding workflow:
    annot_magma_batch.py
    annot_smr_batch.py
Requires following inputs:
    MAGMA GSA outputs
    SMR summary stats
'''

from gettext import find


def main(args):
    from fnmatch import fnmatch
    force = '-f' if args.force else ''
    
    # array submitter
    from ._utils.slurm import array_submitter
    submitter = array_submitter(
        name = f'annot_manhattan_{args.pheno[0]}',
        timeout = 60,
        # debug = True
    )
    
    from ._utils.path import find_gwas
    pheno = find_gwas(args.pheno, dirname = args.gwa, ext = 'fastGWA', long = True)
    for g, p in pheno:
        # make output directory
        if not os.path.isdir(f'{args.out}/{g}'): os.system(f'mkdir -p {args.out}/{g}')
        submitter.add(
            f'python annot_manhattan.py {g} --prefix {p} ' +
            f'--magma {args.magma} --smr {args.smr} --ref {args.ref} -o {args.out} {force}'
        )
    submitter.submit()
        
if __name__ == '__main__':
    from ._utils.slurm import slurm_parser
    parser = slurm_parser(description = 'Plots Manhattan plots for gene-level statistics')
    parser.add_argument('pheno', help = 'Phenotype group', nargs = '*')
    parser.add_argument('--gwa', help = 'Directory to GWAS sumstats, for directory scanning',
      default = '../gwa')
    parser.add_argument('--magma', help = 'MAGMA output directory',
      default = '../annot/magma')
    parser.add_argument('--smr', help = 'SMR output directory',
      default = '../annot/smr')
    parser.add_argument('-r', '--ref', dest = 'ref', help = 'Gene label document',
      default = '../params/genes_ref.txt')
    parser.add_argument('-o','--out', dest = 'out', help = 'Output directory',
      default = '../annot/manhattan/')
    parser.add_argument('-f','--force',dest = 'force', help = 'force overwrite',
      default = False, action = 'store_true')
    args = parser.parse_args()
    # path normalisation
    args.pheno.sort()
    import os
    for arg in ['gwa','magma','smr','ref','out']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
        
    from ._utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args.smr, __file__)
    proj.add_output(args.magma, __file__)
    proj.add_output(args.out, __file__)
    try: main(args)
    except: cmdhistory.errlog()