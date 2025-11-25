'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2025-11-17

A wrapper script to conduct GNOVA correlation analysis

Requires following inputs:
    LDSC-munged GWAS summary statistics
    GNOVA-formatted gene sets and ldscores
Outputs:
    GNOVA report
'''

import os, warnings

def main(args):
    from ._utils.slurm import array_submitter
    submitter = array_submitter(name = 'gcorr_gnova_'+'_'.join(args.p1)+'_'+'_'.join(args.p2),
        n_cpu = 4, timeout = 30, env = args.gnova)
    
    # find and pair GWAS files
    from ._utils.path import find_gwas, pair_gwas
    exposures = find_gwas(args.p1, dirname = args._in, ext = 'sumstats', long = True)
    outcomes = find_gwas(args.p2, dirname = args._in, ext = 'sumstats', long = True)
    pairwise = pair_gwas(exposures, outcomes)

    # find annot files
    gene_sets = []
    for x in os.listdir(args.gene_set):
        if x.startswith('.') or not os.path.isdir(f'{args.gene_set}/{x}'): continue
        if not any([f'{chrom}.gnova' in os.listdir(f'{args.gene_set}/{x}') for chrom in range(1,23)]): 
            warnings.warn(f'Missing GNOVA annotation file in {args.gene_set}/{x}, skipping')
            continue
        gene_sets.append(x)

    print(f'Found following gene sets for GNOVA analysis: \n' + '\n'.join(gene_sets))

    for g1, p1, g2, p2 in pairwise:
        if g1 > g2 or (g1 == g2 and p1 > p2):
            g1, p1, g2, p2 = g2, p2, g1, p1
        outdir = f'{args.out}/{g1}.{g2}'
        os.makedirs(outdir, exist_ok = True)
        for gene_set in gene_sets:
            out_prefix = f'{outdir}/{g1}_{p1}.{g2}_{p2}.{gene_set}.gnova'
            ld_file = f'{args.gene_set}/{gene_set}.gnova.ldscore'
            if (not args.force) and os.path.isfile(f'{out_prefix}'): continue
            cmd = ['python',f'{args.gnova}/GNOVA/gnova.py', '--bfile', args.bfile,
                    '--annot', f'{args.gene_set}/{gene_set}/@.gnova',
                    '--use-ld' if os.path.isfile(f'{ld_file}.csv.gz') else '--save-ld', ld_file,
                    '--out', out_prefix,
                    f'{args._in}/{g1}/{p1}.sumstats', f'{args._in}/{g2}/{p2}.sumstats'
                   ]
            submitter.add(' '.join(cmd))
    submitter.submit()
    return submitter

if __name__ == '__main__':
    from ._utils.slurm import slurm_parser
    parser = slurm_parser(description = 'A wrapper script to conduct GNOVA correlation analysis')
    parser.add_argument('-p1', '--p1', nargs = '+', help = 'Phenotype group 1')
    parser.add_argument('-p2', '--p2', nargs = '*', default = [], help = 'Phenotype group 2 (if not specified, estimate pairwise correlation of p1)')
    parser.add_argument('-i','--in', dest = '_in', help = 'Input directory containing LDSC-munged GWAS summary statistics',
        default = '../gcorr/ldsc_sumstats/')
    parser.add_argument('-o', '--out', help = 'Output directory',
        default = '../gcorr/gnova/')
    parser.add_argument('--gnova', help = 'GNOVA executable directory',
        default = '/rds/project/rds-Nl99R8pHODQ/toolbox/gnova')
    parser.add_argument('--gene_set', help = 'Path to GNOVA-formatted gene set file',
        default = '/rds/project/rds-Nl99R8pHODQ/multiomics/snp_annot')
    parser.add_argument('--bfile', help = 'Path to PLINK binary files for reference panel',
        default = '/rds/project/rds-Nl99R8pHODQ/ref/1000g_eur_ldsc/chr@')
    parser.add_argument('-f','--force', action = 'store_true', help = 'Force overwrite')

    args = parser.parse_args()
    import os
    for arg in ['_in','out','gnova','gene_set']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    
    from ._utils import cmdhistory, logger
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()
