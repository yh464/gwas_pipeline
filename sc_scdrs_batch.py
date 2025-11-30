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

def generate_gene_sets(gene_sumstats, args):
    # generate gene sets from the gene-level summary stats
    import pandas as pd
    from scipy.stats import false_discovery_control as fdr

    df = pd.read_table(gene_sumstats, sep = '\\s+').sort_values(args.pcol, ascending = True).reset_index(drop = True)
    n_genes = df.shape[0]
    nmin = args.nmin * n_genes if 0 < args.nmin < 1 else args.nmin
    nmin = int(max(50, nmin))
    nmax = args.nmax * n_genes if 0 < args.nmax < 1 else args.nmax
    nmax = int(min(nmax, 0.2 * n_genes))
    if nmax <= nmin: raise ValueError('Too few genes in dataset')
    if args.pval == None:
        q = fdr(df[args.pcol])
        nsig = sum(q < 0.05)
    else:
        nsig = sum(df[args.pcol] < args.pval)
    # nsig = max(nsig, nmin); nsig = min(nsig, nmax)
    gene_list = df[args.gcol].tolist()
    if args.bcol != None: gene_weight = df[args.bcol].values
    else: gene_weight = 1
    return pd.DataFrame(dict(gene = gene_list, weight = gene_weight)), [100, 200, 500, 1000, 2000, nsig]

def main(args):
    import os
    import pandas as pd
    from _utils.path import find_gwas, find_gene_sumstats
    from _utils.slurm import array_submitter
    pheno = find_gwas(args.pheno, long = True)
    submitter = array_submitter(name = 'sc_scdrs_'+'_'.join(args.pheno), n_cpu = 4, timeout = 720, env = 'gentoolspy')

    # scans hdf5 files
    for g, p in pheno:
        os.makedirs(f'{args.out}/{g}/{p}', exist_ok = True)
        # check annotation file exists
        annot = find_gene_sumstats(g, p, args._in, args.annot)
        if not annot: Warning(f'Missing gene-level sumstats for {g}/{p}'); continue
        # generate gene set
        weights_file = f'{args.out}/{g}/{p}.genes.txt'
        weights, nsigs = generate_gene_sets(annot, args)
        weights.to_csv(weights_file, index = False, sep = '\t')
        
        for sc in args.sc:
            h5ad_dir = f'{args.h5ad}/{sc}'
            if not os.path.isdir(h5ad_dir): continue
            outdir = f'{args.out}/{g}/{p}/{sc}'
            os.makedirs(outdir, exist_ok = True)
            h5ad = [f'{h5ad_dir}/{x}' for x in os.listdir(h5ad_dir) if x[-5:] =='.h5ad']
            h5ad_prefix = [x[:-5] for x in os.listdir(h5ad_dir) if x[-5:] =='.h5ad']

            for h5, h5prefix in zip(h5ad,h5ad_prefix):
                out_prefix = f'{outdir}/{p}.{h5prefix}.scdrs'
                if not os.path.isfile(f'{out_prefix}.sensitivity.txt'):
                    submitter.add(
                        f'python sc_scdrs_sensitivity.py -i {weights_file} --h5ad {h5} --label '+
                        ' '.join(args.label)+
                        f' -n {" ".join([str(x) for x in nsigs])} -o {out_prefix}'+
                        (' -f' if args.force else '')
                    )
                if os.path.isfile(f'{out_prefix}.score.txt') and os.path.isfile(f'{out_prefix}.enrichment.txt') \
                    and os.path.isfile(f'{out_prefix}.score.png') and \
                    (not args.downstream or (os.path.isfile(f'{out_prefix}.downstream.txt') and \
                        os.path.isfile(f'{out_prefix}.downstream.revigo.txt') and os.path.isfile(f'{args.out}.pseudotime.png'))) \
                    and not args.force: continue
                cmd = ['python', 'sc_scdrs.py', '-i', weights_file, '-n', 
                       f'{nsigs[-1]}' if args.nsig < 0 else f'{args.nsig:.0f}',
                       '--h5ad', h5, '--label'] + args.label
                if args.pval != None: cmd.append(f'-p {args.pval}')
                if args.downstream: cmd.append('-d')
                cmd += ['-o', out_prefix]
                if args.force: cmd.append('-f')
                submitter.add(' '.join(cmd))    
    submitter.submit()
    return submitter

if __name__ == '__main__':  
    from _utils.slurm import slurm_parser
    parser = slurm_parser(description = 'This script runs cell-type enrichments using scDRS')
    parser.add_argument('pheno', nargs = '*', help = 'Phenotypes')
    parser.add_argument('-s','--sc', nargs = '*', help = 'single-cell dataset', default = ['siletti_2023','wang_2025','keefe_2025'])
    parser.add_argument('-i','--in', dest = '_in', help = 'Directory containing gene-level summary statistics',
        default = '../annot/magma')
    parser.add_argument('--annot', help = 'Annotation used to generate gene-level sumstats', default = 'ENSG')
    parser.add_argument('--pcol', help = 'P-value column of input summary stats', default = 'P') # p_SMR
    parser.add_argument('--bcol', help = 'Effect size column of input summary stats', default = 'ZSTAT') # b_SMR
    parser.add_argument('--gcol', help = 'Gene ID column of input summary stats', default = 'GENE') # probeID
    parser.add_argument('-p','--pval', help = 'P-value threshold for significant genes, default FDR 0.05', default = None)
    parser.add_argument('--nmin', help = 'Minimum number of significant genes, number or fraction', type = float, default = 100)
    parser.add_argument('--nmax', help = 'Maximum number of significant genes, number or fraction', type = float, default = 0.1)
    parser.add_argument('-n','--nsig', help = 'Number of significant genes, -1 for all FDR-significant genes, default 1000', 
        type = int, default = 1000)
    parser.add_argument('--h5ad', help = 'Input directory containing h5ad single-cell multiomics dataset',
        default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/multiomics/scdrs') # intentionally absolute
    parser.add_argument('--label', nargs = '*', help = 'Columns containing cell classifications/annotations in the h5ad dataset',
        default = ['ROIGroup', 'ROIGroupCoarse', 'ROIGroupFine', 'roi', 'supercluster_term', 'cluster_id', 'subcluster_id', 'development_stage', # siletti
        'Class','Subclass','Type_updated', 'Cluster', 'Tissue', # wang
        'subcluster_identity_broad','subcluster_identity', # keefe
        ])
    parser.add_argument('-d', '--downstream', help = 'Conduct downstream analyses', default = False, action = 'store_true')
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