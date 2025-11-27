#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2025-11-27

Wrapper script to submit jobs of mix3r

Requires following inputs: 
    GWAS summary statistics, ldsc-munged
'''

import os
from _utils.path import find_gwas
from _utils.gadgets import namespace
from _utils.slurm import array_submitter

def main(args = None, **kwargs):
    if args == None: args = namespace(**kwargs)

    pheno = find_gwas(args.pheno, dirname = args._in, long = True, ext = 'sumstats')
    if len(pheno) != 3:
        raise ValueError(f'Expecting 3 GWAS files for mix3r, found {len(pheno)}')
    os.makedirs(args.out, exist_ok = True)
    out_file = f'{args.out}/{pheno[0][0]}_{pheno[0][1]}.{pheno[1][0]}_{pheno[1][1]}.{pheno[2][0]}_{pheno[2][1]}.json'

    if os.path.isfile(out_file) and not args.force: return

    tmpdir = '/rds/user/yh464/hpc-work/tmp/mix3r'
    os.makedirs(tmpdir, exist_ok = True)
    config = f'{tmpdir}/{pheno[0][0]}_{pheno[0][1]}.{pheno[1][0]}_{pheno[1][1]}.{pheno[2][0]}_{pheno[2][1]}.config.json'
    with open(config, 'w') as f:
        print(f'''
{{
    "sumstats": [
        "{args._in}/{pheno[0][0]}/{pheno[0][1]}.sumstats",
        "{args._in}/{pheno[1][0]}/{pheno[1][1]}.sumstats",
        "{args._in}/{pheno[2][0]}/{pheno[2][1]}.sumstats"
    ],
    "template_dir": "{args.mix3r}/template",
    "nbin_het_hist": 64,

    "out": "{out_file}",

    "snp_filters": {{
        "chromosomes": [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22],
        "maf_thresh": 0.05,
        "info_thresh": 0.8,
        "z_thresh": 32,
        "exclude_regions": ["6:25000000-34000000"]
    }},

    "pruning": {{
        "do_pruning": true,
        "r2_prune_thresh": 0.8,
        "n_random": 300000,
        "rand_prune_seed": 1
    }},

    "optimization": {{
        "maxiter_1d_glob": 128,
        "maxiter_1d_loc": 200,
        "maxiter_2d_glob": 128,
        "maxiter_2d_loc": 200,
        "maxiter_3d": 16
    }}
}}
''', file = f)
    
    cmd = f'python {args.mix3r}/mix3r_int_weights.py --config {config}'
    submitter = array_submitter('mix3r', partition = 'sapphire', n_cpu = 6, timeout = 720,
        env = args.mix3r)
    submitter.add(cmd)
    submitter.submit()

if __name__ == '__main__':
    from _utils.slurm import slurm_parser
    parser = slurm_parser(description = 'This script runs mix3r on 3 GWAS summary statistics')
    parser.add_argument('pheno', type = str, nargs = '*', help = 'Phenotypes, exactly 3 traits')
    parser.add_argument('-i', '--in', dest = '_in', type = str, help = 'Input directory',
        default = '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/gcorr/ldsc_sumstats')
    parser.add_argument('--mix3r', type = str, help = 'Mix3r installation path', 
        default = '/rds/project/rds-Nl99R8pHODQ/toolbox/mix3r')
    parser.add_argument('-o', '--out', type = str, help = 'Output directory',
        default = '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/gcorr/mix3r')
    parser.add_argument('-f', '--force', action = 'store_true', help = 'Force overwrite')
    args = parser.parse_args()
    
    for attr in ['_in', 'out', 'mix3r']:
        setattr(args, attr, os.path.abspath(getattr(args, attr)))

    from _utils import cmdhistory, logger
    logger.splash(args); cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()
    