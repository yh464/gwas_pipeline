#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2024-11-12
Version 2: 2025-08-22

Conducts MAGMA GSEA for all phenotypes

Requires following inputs: 
    MAGMA annotation files (including H-MAGMA)
    MAGMA binary
'''
def main(args = None, **kwargs):
    from _utils.gadgets import namespace
    if args == None:
        from _utils.gadgets import namespace
        args = namespace(**kwargs)

    from _utils.slurm import array_submitter
    submitter = array_submitter('sc_magma_gsea_'+'_'.join(args.pheno),
        n_cpu = 1, timeout = 20, parallel = 4, modules = ['gcc/11'], env = 'base')

    # find gene set files
    import os
    gsets = [(f'{args.gset}/{x}', x[:-4]) for x in os.listdir(args.gset) if x[-4:] == '.txt']
    gscores = [(f'{args.gscore}/{x}', x[:-4]) for x in os.listdir(args.gscore) if x[-4:] == '.txt']

    # identify phenotypes
    from _utils.path import find_gwas, find_gene_sumstats
    pheno = find_gwas(args.pheno, long = True)

    for g, p in pheno:
        outdir = f'{args.out}/{g}/{p}'
        os.makedirs(outdir, exist_ok = True)

        sumstat = find_gene_sumstats(g, p, args._in, args.annot, ext = 'genes.raw')
        if not sumstat: Warning(f'No gene-level sumstat found for {g}/{p}'); continue

        for gset, gset_prefix in gsets:
            out_prefix = f'{outdir}/{p}.{args.annot}.{gset_prefix}'
            if os.path.isfile(f'{out_prefix}.gsa.out') and not args.force: continue
            submitter.add(f'{args.magma} --gene-results {sumstat} --set-annot {gset} --out {out_prefix}')
        
        for gscore, gscore_prefix in gscores:
            out_prefix = f'{outdir}/{p}.{args.annot}.{gscore_prefix}'
            if os.path.isfile(f'{out_prefix}.gsa.out') and not args.force: continue
            submitter.add(f'{args.magma} --gene-results {sumstat} --gene-covar {gscore} --out {out_prefix}')
    submitter.submit()
    return submitter

if __name__ == '__main__':
    from _utils.slurm import slurm_parser
    parser = slurm_parser(description = 'This script runs cell-type/pathway enrichments using MAGMA GSEA')
    parser.add_argument('pheno', nargs = '*', help = 'Phenotypes')
    parser.add_argument('-i','--in', dest = '_in', help = 'Directory containing gene-level summary statistics',
        default = '../annot/magma')
    parser.add_argument('--annot', help = 'Annotation used to generate gene-level sumstats', default = 'ENSG')
    parser.add_argument('--gset', dest = 'gset', help = 'Gene sets to study enrichment, scans directory',
        default = '../multiomics/gene_set')
    parser.add_argument('--gscore', help = 'Directory containing gene scores', default = '../multiomics/gene_score')
    parser.add_argument('--magma', dest = 'magma', help = 'MAGMA executable', default = '../toolbox/magma/magma')
    parser.add_argument('-o', '--out', dest = 'out', help = 'output directory', default = '../sc/magma_gsea')
    parser.add_argument('-f','--force',dest = 'force', help = 'force overwrite', default = False, action = 'store_true')
    args = parser.parse_args()
    
    # path normalisation
    import os
    for arg in ['_in','gset','gscore','magma','out']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))

    from _utils import cmdhistory, logger
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()