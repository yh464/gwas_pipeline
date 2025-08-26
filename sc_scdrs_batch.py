#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2025-08-21

Conducts scDRS scoring for a group of phenotypes in a job batch

Requires following inputs: 
    Gene-level summary statistics
        Column specifications: gene ID (ENSGxxxxxxxx), p-value and effect size (optional)
    P-value threshold for sig genes (default fdr < 0.05)
    Minimum and maximum number of genes to include as sig genes (default: min 100, max 10% of genes)
    Pre-processed scDRS single cell data
        Column specifications: Columns of adata.obs that contain cell group labels
            (for Siletti et al 2023: ROIGroup, ROIGroupCoarse, ROIGroupFine, roi, supercluster_term, cluster_id, subcluster_id, development_stage)
'''

def main(args):
    import os
    import pandas as pd
    from _utils.path import find_gwas, find_gene_sumstats
    from _utils.slurm import array_submitter
    pheno = find_gwas(args.pheno, long = True)
    submitter = array_submitter(name = 'sc_scdrs_'+'_'.join(args.pheno), n_cpu = 16, timeout = 360, env = 'gentoolspy')

    # scans hdf5 files
    h5ad = [f'{args.h5ad}/{x}' for x in os.listdir(args.h5ad) if x[-5:] =='.h5ad']
    h5ad_prefix = [x[:-5] for x in os.listdir(args.h5ad) if x[-5:] =='.h5ad']

    for g, p in pheno:
        outdir = f'{args.out}/{g}/{p}'
        if not os.path.isdir(outdir): _ = os.system(f'mkdir -p {outdir}')

        for h5, h5prefix in zip(h5ad,h5ad_prefix):
            # check annotation file exists
            annot = find_gene_sumstats(g,p, args._in, args.annot)
            if not annot: Warning(f'Missing gene-level sumstats for {g}/{p}'); continue

            out_prefix = f'{outdir}/{p}.{h5prefix}.scdrs'
            if os.path.isfile(f'{out_prefix}.score.txt') and os.path.isfile(f'{out_prefix}.enrichment.txt') \
                and not args.force: continue
            cmd = ['python', 'sc_scdrs.py', '-i', annot, '--pcol', args.pcol, '--bcol', args.bcol,
                   '--gcol', args.gcol, f'--nmax {args.nmax} --nmin {args.nmin}', '--h5ad', h5, '--label'] + args.label
            if args.pval != None: cmd.append(f'-p {args.pval}')
            cmd += ['-o', out_prefix]
            if args.force: cmd.append('-f')
            submitter.add(' '.join(cmd))
        if all([os.path.isfile(f'{outdir}/{p}.{h5prefix}.scdrs.enrichment.txt') for h5prefix in h5ad_prefix]):
            df = pd.concat([pd.read_table(f'{outdir}/{p}.{h5prefix}.scdrs.enrichment.txt').assign(dataset = h5prefix)
                for h5prefix in h5ad_prefix], axis = 0)
            df.to_csv(f'{args.out}/{g}/{p}.scdrs.enrichment.txt', index = False, sep = '\t')
    submitter.submit()
    return submitter

if __name__ == '__main__':  
    from _utils.slurm import slurm_parser
    parser = slurm_parser(description = 'This script runs cell-type enrichments using scDRS')
    parser.add_argument('pheno', nargs = '*', help = 'Phenotypes')
    parser.add_argument('-i','--in', dest = '_in', help = 'Directory containing gene-level summary statistics',
        default = '../annot/magma')
    parser.add_argument('--annot', help = 'Annotation used to generate gene-level sumstats', default = 'ENSG')
    parser.add_argument('--pcol', help = 'P-value column of input summary stats', default = 'P') # p_SMR
    parser.add_argument('--bcol', help = 'Effect size column of input summary stats', default = 'ZSTAT') # b_SMR
    parser.add_argument('--gcol', help = 'Gene ID column of input summary stats', default = 'GENE') # probeID
    parser.add_argument('-p','--pval', help = 'P-value threshold for significant genes, default FDR 0.05', default = None)
    parser.add_argument('--nmin', help = 'Minimum number of significant genes, number or fraction', type = float, default = 100)
    parser.add_argument('--nmax', help = 'Maximum number of significant genes, number or fraction', type = float, default = 0.1)
    parser.add_argument('--h5ad', help = 'Input directory containing h5ad single-cell multiomics dataset',
        default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/multiomics/scdrs/siletti_2023') # intentionally absolute
    parser.add_argument('--label', nargs = '*', help = 'Columns containing cell classifications/types in the h5ad dataset',
        default = ['ROIGroup', 'ROIGroupCoarse', 'ROIGroupFine', 'roi', 'supercluster_term', 'cluster_id', 'subcluster_id', 'development_stage'])
    parser.add_argument('-o', '--out', dest = 'out', help = 'output directory', default = '../sc/scdrs')
    parser.add_argument('-f','--force',dest = 'force', help = 'force overwrite', default = False, action = 'store_true')
    args = parser.parse_args()
    
    # path normalisation
    import os
    for arg in ['_in','h5ad','out']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))

    from _utils import logger, cmdhistory
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()