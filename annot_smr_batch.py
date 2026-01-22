#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-01-16
Version 2: 2026-01-22

Conducts Summary-data Mendelian Randomisation to prioritise genes using
expression, isoform, splicing, methylation and chromosome accessibility QTL.

Requires following inputs: 
    GWAS summary statistics,
    BESD format xQTL datasets.
'''

from _utils.logger import logger
import pandas as pd
import numpy as np
import os
from fnmatch import fnmatch
log = logger()

def format_gwa(gwa, tmpgwa):
    # input file name, usually fastGWA format
    with open(gwa) as f:
        hdr = f.readline().strip().split()
    out = [None, None, None, None, None, None, None, None] # 8 columns
    out_col = [None] * 8
    os.makedirs(os.path.dirname(tmpgwa), exist_ok = True)
    
    for idx, col in enumerate(hdr):
        if col.lower() in ['snp','id','rsid']: out[0] = idx; out_col[0] = col
        if col.lower() in ['a1','ref','refallele','effectallele']: out[1] = idx; out_col[1] = col
        if col.lower() in ['a2','alt','altallele','otherallele']: out[2] = idx; out_col[2] = col
        if col.lower() in ['af1','freq','eaf','maf']: out[3] = idx; out_col[3] = col
        if col.lower() in ['beta','logor','b']: out[4] = idx; out_col[4] = col; log_eff = False
        if col.lower() in ['or']: out[4] = idx; out_col[4] = col; log_eff = True
        if col.lower() in ['se', 'stderr']: out[5] = idx; out_col[5] = col
        if col.lower() in ['p','pval']: out[6] = idx; out_col[6] = col
        if col.lower() in ['nobs','n']: out[7] = idx; out_col[7] = col
    
    log.log(f'''Identified necessary columns at positions: 
            SNP  = {out[0]}
            A1   = {out[1]}
            A2   = {out[2]}
            AF1  = {out[3]}
            BETA = {out[4]}
            SE   = {out[5]}
            P    = {out[6]}
            N    = {out[7]}''')

    if any([x is None for x in out[:3]]): raise ValueError('Missing necessary columns')

    if not log_eff:
        print_field = ','.join([f'${x+1}' for x in out])
        cmd = ['awk', '-v', r'OFS="\t"', '\'{print', print_field+'}\'', gwa, '>', tmpgwa]
        log.log(' '.join(cmd))
        os.system(' '.join(cmd))
    else:
        df = pd.read_table(gwa, usecols = out)
        df = df.loc[:, out_col]
        df.iloc[:,4] = np.log(df.iloc[:,4])
        df.to_csv(tmpgwa, sep = '\t', index = False)

def main(args):
    import os
    from fnmatch import fnmatch
    
    # temp directory
    tmpdir = '/home/yh464/rds/rds-rb643-ukbiobank2/Data_Users/yh464/temp/smr_temp'
    if not os.path.isdir(tmpdir): os.system(f'mkdir -p {tmpdir}')
    
    # array submitter
    from _utils.slurm import array_submitter
    submitter = array_submitter(name = f'annot_smr_{args.pheno[0]}',n_cpu = 2,timeout = 90)
    
    # find QTL files
    log.log('Following QTL have been found:')
    qtl_list = []
    for y in os.listdir(args.qtl):
        if fnmatch(y, '*.besd'):
            log.log(y)
            qtl_list.append([f'{args.qtl}/{y}'.replace('.besd','')] * 24)
        if os.path.isdir(f'{args.qtl}/{y}'):
            if any([fnmatch(z, '*.besd') for z in os.listdir(f'{args.qtl}/{y}')]): log.log(y)
            else: continue
            tmp_list = []
            for chrom in range(1, 25):
                found = False
                for f in os.listdir(f'{args.qtl}/{y}'):
                    if fnmatch(f.replace('X','23').replace('Y','24'), f'*chr{chrom}.besd'):
                        tmp_list.append(f'{args.qtl}/{y}/'+f.replace('.besd',''))
                        found = True
                        break
                if found: continue
            qtl_list.append(tmp_list)

    # parse input PLINK binaries
    if os.path.isfile(f'{args.bfile}.bed'):
        bfile_list = [args.bfile] * 24
    else:
        bfile_list = []
        for chrom in range(1,25):
            found = False
            for f in os.listdir(args.bfile):
                if fnmatch(f.replace('X','23').replace('Y','24'), f'*chr{chrom}.bed'):
                    bfile_list.append(f'{args.bfile}/'+f.replace('.bed',''))
                    found = True
                    break
            if found: continue
            # if there is no match
            bfile_list.append(None)
    
    from _utils.path import find_gwas
    pheno = find_gwas(args.pheno, dirname = args._in, ext = 'fastGWA', long = True)
    for g,p in pheno:
        os.makedirs(f'{args.out}/{g}', exist_ok = True)
        tmpgwa = f'{tmpdir}/{g}/{p}.txt'
        gwa = f'{args._in}/{g}/{p}.fastGWA'
        if not os.path.isfile(tmpgwa) or args.force:
            try: format_gwa(gwa, tmpgwa)
            except: log.warn(f'{p} missing necessary columns'); continue

        for qtl in qtl_list:
            if qtl[0] == qtl[1] and bfile_list[0] == bfile_list[1]:
                if os.path.isfile(f'{args.out}/{p}.smr') and not args.force: continue
                submitter.add(f'{args.smr} --bfile {bfile_list[0]} --gwas-summary {tmpgwa} '+
                    f'--beqtl-summary {qtl[0]} --out {args.out}/{p}.{os.path.basename(qtl[0])}')
            else:
                os.makedirs(f'{args.out}/{p}.{os.path.basename(qtl[0])}', exist_ok = True)
                for q, b, chrom in zip(qtl, bfile_list, range(1,25)):
                    if os.path.isfile(f'{args.out}/{p}.{os.path.basename(qtl[0])}/chr{chrom}.smr') and not args.force:
                        continue
                    submitter.add(f'{args.smr} --bfile {b} --gwas-summary {tmpgwa} '+
                        f'--beqtl-summary {q} --out {args.out}/{p}.{os.path.basename(qtl[0])}/chr{chrom}')
    
    submitter.submit()

if __name__ == '__main__':
    from _utils.slurm import slurm_parser
    parser = slurm_parser(
      description = 'This programme batch runs summary data randomisation')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
    parser.add_argument('-i','--in', dest = '_in', help = 'Directory containing all GWA summary statistics',
      default = '../gwa/')
    parser.add_argument('-q','--qtl', dest = 'qtl', help = 'Directory containing all xQTL files',
      default = '../params/xqtl')
    parser.add_argument('-s','--smr', dest = 'smr', help = 'Location of SMR binary',
      default = '/home/yh464/rds/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/smr') # intentionally absolute
    parser.add_argument('-b', '--bfile', dest = 'bfile', help = 'bed binary to use in magma',
      default = '/rds/project/rds-Nl99R8pHODQ/UKB/Imaging_genetics/yh464/bed/') # intentionally absolute
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
      default = '../annot/smr')
    parser.add_argument('-f','--force',dest = 'force', help = 'force overwrite',
      default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','qtl','smr','bfile']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    
    from _utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheng/%pheno_%maf.fastGWA', __file__)
    proj.add_output(args.out+'/%pheng/%pheno.*.smr',__file__)
    try: main(args)
    except: cmdhistory.errlog()