#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2026-02-09

Python wrapper to run pathway-specific PRS using PRSet#

Inputs:
    GWAS summary statistics in fastGWA format (columns: SNP, CHR, POS, A1, A2, BETA/OR, P)
    PLINK .bed files of target population
    Ensembl GTF reference file (defaults to GRCh37)
    MSigDB files containing pathway information (defaults to c3.all = regulatory targets, and c5.go = gene ontology)
'''

import os
from _utils.path import find_gwas, find_bed
from _utils import logger, cmdhistory
from _utils.slurm import array_submitter, slurm_parser
log = logger.logger() 

def main(args):
    pheno = find_gwas(args.pheno, dirname = args._in, long = True)
    bed = find_bed(args.bed)[0].replace('chr1','chr#') # wildcard for chromosomes = #
    submitter = array_submitter('prs_pathway_' + '_'.join([p for _,p in pheno]), n_cpu = 4, timeout = 240,
        partition = 'sapphire')
    gmt_prefix = []
    files_selected = []
    files = [x for x in os.listdir(args.msigdb) if x.endswith('.symbols.gmt')]
    msg = '\n'.join([f'{i}: {file}' for i, file in enumerate(files)] + 
        ['', 'Please select the MSigDB files to use for pathway-specific PRS calculation, empty line to finish.','',
         'Default: c3.all.v2026.1.Hs.symbols.gmt, c5.go.v2026.1.Hs.symbols.gmt', '',''] +
        ['Currently selected:',''] + files_selected + ['\n'])
    while len(fid := input(msg)) > 0:
        try: 
            fid = int(fid)
            files_selected.append(files[fid])
            gmt_prefix.append('.'.join(files[fid].split('.')[0:2]).replace('.all',''))
            os.system('clear')
        except: os.system('clear')
    if len(files_selected) == 0:
        files_selected = ['c3.all.v2026.1.Hs.symbols.gmt', 'c5.go.v2026.1.Hs.symbols.gmt']
        gmt_prefix = ['c3', 'c5.go']
    gmt_prefix = '_'.join(gmt_prefix)
    log.log(f'Selected following MSigDB files for pathway-specific PRS calculation:')
    for file in files_selected: log.log(f'    {file}')

    for file in files_selected:
        os.system(f'cat {args.msigdb}/{file} >> {args.out}/to_analyse.gmt')

    for g, p in pheno:
        os.makedirs(f'{args.out}/{g}/{p}', exist_ok = True)
        target_file = f'{args.out}/{g}/{p}/{p}_prset_{gmt_prefix}.best'
        if not args.force and os.path.isfile(target_file): continue
        cmd = [f'{args.prsice}/bin/PRSice', '--base', f'{args._in}/{g}/{p}.fastGWA', '--target', bed, '--ld', args.ref,
               '--out', f'{args.out}/{g}/{p}/{p}_prset_{gmt_prefix}',
               '--msigdb', f'{args.out}/to_analyse.gmt', '--gtf', args.gtf
               ]
        if 'OR' in open(f'{args._in}/{g}/{p}.fastGWA').readline(): cmd.extend(['--or'])
        submitter.add(' '.join(cmd))
    submitter.submit()

if __name__ == '__main__':
    parser = slurm_parser(description = 'Run pathway-specific PRS using PRSet')
    parser.add_argument('pheno', nargs = '*', help = 'Phenotypes to analyse pathway PRS')
    parser.add_argument('-i','--in', dest = '_in', default = '../gwa', help = 'Input directory containing GWAS summary statistics in fastGWA format')
    parser.add_argument('-o','--out', default = '../prs', help = 'Output directory')
    parser.add_argument('--prsice', default = '/rds/project/rds-Nl99R8pHODQ/toolbox/PRSice', 
        help = 'Path to PRSice executable directory') # intentionally absolute
    parser.add_argument('--bed', default = '/rds/project/rds-Nl99R8pHODQ/UKB/Imaging_genetics/yh464/bed/', 
        help = 'Path to PLINK .bed file of target population') # intentionally absolute
    parser.add_argument('--ref', default = '/rds/project/rds-Nl99R8pHODQ/UKB/Imaging_genetics/yh464/bed/chr#',
        help = 'Reference LD panel in PLINK format') # intentionally absolute, with wildcard # for chromosome number
    parser.add_argument('--gtf', default = '/rds/project/rds-Nl99R8pHODQ/ref/ensg/ensg.*build*.gtf.txt',
        help = 'Path to GTF reference file') # intentionally absolute, with build specified in file name
    parser.add_argument('--build', default = 'hg19', choices = ['hg19', 'hg38'], 
        help = 'Genome build of the reference GTF file, default to hg19')
    parser.add_argument('--msigdb', default = '/rds/project/rds-Nl99R8pHODQ/ref/msigdb/', 
        help = 'Path to MSigDB files, default to /rds/project/rds-Nl99R8pHODQ/ref/msigdb/')
    parser.add_argument('-f','--force', action = 'store_true', help = 'Whether to overwrite existing output files')
    args = parser.parse_args()

    args.gtf = args.gtf.replace('*build*', args.build)
    for arg in ['_in', 'out', 'prsice', 'bed', 'ref', 'gtf', 'msigdb']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    
    cmdhistory.log()
    logger.splash(args)
    try: main(args)
    except: cmdhistory.errlog()
