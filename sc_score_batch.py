#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2025-08-28

Generates gene scores for MAGMA gene covariate analysis

Requires following inputs: 
    directory containing single-cell h5ad files
    columns containing cell classifications
'''

def main(args):
    from ._utils.slurm import array_submitter
    submitter = array_submitter(name = 'sc_score', n_cpu = 16, timeout = 480, env = 'gentoolspy')

    # finds h5ad files
    import os
    os.makedirs(args.out, exist_ok = True)
    h5ad_files = [f'{args._in}/{x}' for x in os.listdir(args._in) if x[-5:] =='.h5ad']
    h5ad_prefix = [x[:-5] for x in os.listdir(args._in) if x[-5:] =='.h5ad']
    if len(h5ad_files) == 0: raise ValueError('No h5ad files found in input directory')

    for h5, prefix in zip(h5ad_files, h5ad_prefix):
        out_file = f'{args.out}/{prefix}.cepo.txt'
        if os.path.isfile(out_file) and not args.force: continue
        # cmd = ['Rscript', 'sc_score.r','-i', h5, '-o', f'{args.out}/{prefix}', '--label'] + args.label
        cmd = ['python', 'sc_score_cepo.py','-i', h5, '-o', f'{args.out}/{prefix}', '--label'] + args.label
        if args.force: cmd.append('-f')
        submitter.add(' '.join(cmd))
    submitter.submit()
    return submitter

if __name__ == '__main__':
    from ._utils.slurm import slurm_parser
    parser = slurm_parser(description = 'This script generates gene scores for MAGMA gene covariate analysis')
    parser.add_argument('-i','--in', dest = '_in', help = 'Input directory containing h5ad single-cell multiomics dataset',
        default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/multiomics/cepo') # intentionally absolute
    parser.add_argument('-o','--out', help = 'Output directory', 
        default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/multiomics/gene_score') # intentionally absolute
    parser.add_argument('--label', nargs = '*', help = 'Columns containing cell classifications/annotations in the h5ad dataset',
        default = ['ROIGroup', 'ROIGroupCoarse', 'ROIGroupFine', 'roi', 'supercluster_term', 'cluster_id', 'subcluster_id', 'development_stage', # siletti
        'Class','Subclass','Type_updated', 'Cluster', 'Tissue']) # wang
    parser.add_argument('-f','--force', help = 'Force overwrite', action = 'store_true', default = False)
    args = parser.parse_args()
    
    import os
    args._in = os.path.realpath(args._in); args.out = os.path.realpath(args.out)

    from ._utils import logger, cmdhistory
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()