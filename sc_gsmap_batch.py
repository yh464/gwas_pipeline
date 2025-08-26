#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2025-08-26

Conducts gsMap spatial LDSC analysis for a group of phenotypes in a job batch

Requires following inputs: 
    LDSC-munged summary statistics
    gsMap processed spatial transcriptomics data folder (need to be preprocessed - see genetics_qc/ref_processing/preprocess_sc_4gsmap.sh)

'''
from math import e


def main(args):
    from _utils.gadgets import mv_symlink
    from _utils.slurm import array_submitter
    submitter = array_submitter('gsmap_'+'_'.join(args.pheno), timeout = 240, n_cpu = 16, env = args.gsmap)

    # find ST datasets
    import os
    st_datasets = os.listdir(args.st)

    from _utils.path import find_gwas
    pheno = find_gwas(args.pheno, dirname = args._in, ext = 'sumstats', long = True)

    for g, p in pheno:
        outdir = f'{args.out}/{g}/{p}'
        os.makedirs(outdir, exist_ok = True)

        for s in st_datasets:
            gsmap_out_ldsc = f'{args.st}/{s}/spatial_ldsc/{s}_{g}_{p}.csv.gz'
            gsmap_out_cauchy = f'{args.st}/{s}/spatial_ldsc/{s}_{g}_{p}.Cauchy.csv.gz'
            gsmap_out_rpt = f'{args.st}/{s}/report/{g}_{p}/gsMap_plot'
            out_ldsc = f'{outdir}/{s}_spatial_ldsc.csv.gz'
            out_cauchy = f'{outdir}/{s}_spatial_ldsc.cauchy.csv.gz'
            out_rpt = f'{outdir}/{s}_gsmap_report'

            cmds = []

            # spatial ldsc
            if not os.path.isfile(gsmap_out_ldsc) or args.force:
                cmds.append(f'gsmap run_spatial_ldsc --workdir {args.st} --sample_name {s} --trait_name {g}_{p} '+
                    f'--sumstats_file {args._in}/{g}/{p}.sumstats --w_file {args.gsmap}/gsMap_resource/LDSC_resource/weights_hm3_no_hla/weights. '+
                    '--num_processes 16')
            mv_symlink(gsmap_out_ldsc, out_ldsc)

            # Cauchy combination
            if not os.path.isfile(gsmap_out_cauchy) or args.force:
                cmds.append(f'gsmap run_cauchy_combination --workdir {args.st} --sample_name {s} --trait_name {g}_{p} --annotation annotation')
            mv_symlink(gsmap_out_cauchy, out_cauchy)

            # report generation
            if args.report:
                os.makedirs(out_rpt, exist_ok = True)
                os.symlink(out_rpt, gsmap_out_rpt) # symlink the gsmap output directory to output filesystem
                if not os.path.isfile(f'{out_rpt}/{s}_{g}_{p}_gsMap_plot.html') or args.force:
                    cmds.append(f'gsmap run_report --workdir {args.st} --sample_name {s} --trait_name {g}_{p} --annotation annotation'+
                                f'--sumstats_file {args._in}/{g}/{p}.sumstats --top_corr_genes 50')
            submitter.add(*cmds) # ensure that these commands are run in the same file
    submitter.submit()


if __name__ == '__main__':  
    from _utils.slurm import slurm_parser
    parser = slurm_parser(description = 'This script runs cell-type enrichments using scDRS')
    parser.add_argument('pheno', nargs = '*', help = 'Phenotypes')
    parser.add_argument('-i','--in', dest = '_in', help = 'Directory containing LDSC summary statistics',
        default = '../gcorr/ldsc_sumstats')
    parser.add_argument('-s','--st', help = 'Directory containing gsMap processed spatial transcriptomics data',
        default = '/rds/project/rds-Nl99R8pHODQ/multiomics/gsmap') # intentionally absolute
    parser.add_argument('--gsmap', help = 'gsMap package and resources directory',
        default = '/rds/project/rds-Nl99R8pHODQ/toolbox/gsmap') # intentionally absolute
    parser.add_argument('-o', '--out', help = 'output directory', default = '../sc/gsmap')
    parser.add_argument('-r','--report', help = 'generate report', default = False, action = 'store_true')
    parser.add_argument('-f','--force', help = 'force overwrite', default = False, action = 'store_true')
    args = parser.parse_args()
    
    # path normalisation
    import os
    for arg in ['_in','st','gsmap','out']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))

    from _utils import logger, cmdhistory
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()