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

def main(args):
    import pandas as pd
    import os
    from _utils.path import find_gwas
    from _utils.gadgets import mv_symlink
    from _utils.slurm import array_submitter
    ldsc_submitter = array_submitter('gsmap_'+'_'.join(args.pheno), timeout = 120, n_cpu = 4, env = args.gsmap)
    cauchy_submitter = array_submitter('gsmap_cauchy_'+'_'.join(args.pheno), timeout = 10, n_cpu = 4, env = args.gsmap, dependency=ldsc_submitter)
    rpt_submitter = array_submitter('gsmap_rpt_'+'_'.join(args.pheno), timeout = 90, n_cpu = 20, env = args.gsmap, dependency=ldsc_submitter)

    # find ST datasets and phenotype files
    st_datasets = os.listdir(args.st)
    pheno = find_gwas(args.pheno, dirname = args._in, ext = 'sumstats', long = True)
    

    for g, p in pheno:
        outdir = f'{args.out}/{g}/{p}'
        os.makedirs(outdir, exist_ok = True)

        for s in st_datasets:
            gsmap_out_ldsc = f'{args.st}/{s}/spatial_ldsc/{s}_{g}_{p}.csv.gz'
            gsmap_out_cauchy = f'{args.st}/{s}/cauchy_combination/{s}_{g}_{p}.Cauchy.csv.gz'
            if args.report: os.makedirs(f'{args.st}/{s}/report', exist_ok = True)
            gsmap_out_rpt = f'{args.st}/{s}/report/{g}_{p}'
            out_ldsc = f'{outdir}/{s}_spatial_ldsc.csv.gz'
            out_cauchy_ct = f'{outdir}/{s}_spatial_ldsc.cell_type.cauchy.csv.gz'
            out_cauchy_region = f'{outdir}/{s}_spatial_ldsc.region.cauchy.csv.gz'
            out_rpt = f'{outdir}/{s}_gsmap_report'
            print(f'gsMap output to:\n  {gsmap_out_ldsc}\n  {gsmap_out_cauchy}\n  {gsmap_out_rpt}')
            print(f'Output will be moved to:\n  {out_ldsc}\n  {out_cauchy_ct}\n  {out_cauchy_region}\n  {out_rpt}')

            cmds = []

            # spatial ldsc
            if not os.path.isfile(gsmap_out_ldsc) or args.force:
                # check chunked outputs
                n_chunks = max([int(x.replace(f'{s}_chunk','')) for x in os.listdir(f'{args.st}/{s}/generate_ldscore') if x.startswith(f'{s}_chunk')]+[0])
                if n_chunks == 0: raise ValueError(f'No chunked ldscore files found for {s} in {args.st}/{s}/generate_ldscore')
                # calculate 100 chunks at a time
                chunks = [(x, min(x+99, n_chunks)) for x in range(1, n_chunks+1, 100)]
                chunk_files = []
                for start_chunk, end_chunk in chunks:
                    gsmap_out_ldsc_chunk = f'{args.st}/{s}/spatial_ldsc/{s}_{g}_{p}_chunk{start_chunk}-{end_chunk}.csv.gz'
                    chunk_files.append(gsmap_out_ldsc_chunk)
                    if not os.path.isfile(gsmap_out_ldsc_chunk) or args.force:
                        ldsc_submitter.add(f'gsmap run_spatial_ldsc --workdir {args.st} --sample_name {s} --trait_name {g}_{p} '+
                                    f'--sumstats_file {args._in}/{g}/{p}.sumstats --w_file {args.gsmap}/gsMap_resource/LDSC_resource/weights_hm3_no_hla/weights. '+
                                    f'--num_processes 4 --chunk_range {start_chunk} {end_chunk}')
                if all([os.path.isfile(x) for x in chunk_files]) and not os.path.isfile(gsmap_out_ldsc):
                    pd.concat([pd.read_csv(x, index_col = False) for x in chunk_files]).to_csv(gsmap_out_ldsc, index = False)
            mv_symlink(gsmap_out_ldsc, out_ldsc)

            # Cauchy combination
            if not os.path.isfile(out_cauchy_region) or not os.path.isfile(out_cauchy_ct) or args.force:
                cmds.append(f'gsmap run_cauchy_combination --workdir {args.st} --sample_name {s} --trait_name {g}_{p} --annotation annotation')
                cmds.append(f'mv {gsmap_out_cauchy} {out_cauchy_ct}')
                cmds.append(f'gsmap run_cauchy_combination --workdir {args.st} --sample_name {s} --trait_name {g}_{p} --annotation region')
                cmds.append(f'mv {gsmap_out_cauchy} {out_cauchy_region}')
                cmds.append(f'ln -s {out_cauchy_region} {gsmap_out_cauchy}') # create a symlink so that report generation can find the cauchy file
            cauchy_submitter.add(*cmds)

            # report generationp
            if args.report:
                os.makedirs(out_rpt, exist_ok = True)
                if not os.path.islink(gsmap_out_rpt): os.symlink(out_rpt, gsmap_out_rpt) # symlink the gsmap output directory to output filesystem
                if not os.path.isfile(f'{out_rpt}/{s}_{g}_{p}_gsMap_plot.html') or args.force:
                    rpt_submitter.add(f'gsmap run_report --workdir {args.st} --sample_name {s} --trait_name {g}_{p} --annotation region '+
                                f'--sumstats_file {args._in}/{g}/{p}.sumstats --top_corr_genes 50')
    ldsc_submitter.submit()
    cauchy_submitter.submit()
    rpt_submitter.submit()
    return rpt_submitter

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