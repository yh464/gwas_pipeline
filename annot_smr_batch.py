#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-01-16

Conducts Summary-data Mendelian Randomisation to prioritise genes using
expression, isoform, splicing, methylation and chromosome accessibility QTL.

Requires following inputs: 
    GWAS summary statistics,
    BESD format xQTL datasets.
'''

def format_gwa(gwa, tmpgwa):
    # input file name, usually fastGWA format
    import pandas as pd
    import numpy as np
    df = pd.read_table(gwa)
    out = [0,0,0,0,0,0,0,0] # 8 columns
    
    for col in df.columns:
        if col.lower() in ['snp','id','rsid']: out[0] = df[col]
        if col.lower() in ['a1','ref','refallele','effectallele']: out[1] = df[col]
        if col.lower() in ['a2','alt','altallele','otherallele']: out[2] = df[col]
        if col.lower() in ['af1','freq','eaf','maf']: out[3] = df[col]
        if col.lower() in ['beta','logor','b']: out[4] = df[col]
        if col.lower() in ['or']: out[4] = np.log(df[col])
        if col.lower() in ['se', 'stderr']: out[5] = df[col]
        if col.lower() in ['p','pval']: out[6] = df[col]
        if col.lower() in ['nobs','n']: out[7] = df[col]
    
    out = pd.concat(out, axis = 1)
    out.to_csv(tmpgwa, sep = '\t', index = False)

def main(args):
    import os
    from fnmatch import fnmatch
    
    # temp directory
    tmpdir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/temp/smr_temp'
    if not os.path.isdir(tmpdir): os.system(f'mkdir -p {tmpdir}')
    
    # array submitter
    from _utils.slurm import array_submitter
    submitter = array_submitter(
        name = f'annot_smr_{args.pheno[0]}',
        n_cpu = 2,
        timeout = 90,
        debug = True
    )
    
    # annotation utility for single fastGWA and single xQTL dataset
    def annot_smr(gwa, xqtl, bfile, out, smr, force):
        # output directory
        prefix = os.path.basename(gwa)
        prefix = '.'.join(prefix.split('.')[:-1])
        if not os.path.isdir(out): os.system(f'mkdir -p {out}')
        
        # munge input summary statistics
        if not os.path.isfile(f'{tmpdir}/{prefix}.txt') or force:
            format_gwa(gwa, f'{tmpdir}/{prefix}.txt')
        
        # parse input xqtl file
        from fnmatch import fnmatch
        qtl = os.path.basename(xqtl).replace('.besd','')
        if os.path.isfile(f'{xqtl}.besd'):
            xqtl_list = [xqtl] * 24
        else:
            xqtl_list = []
            for chrom in range(1,25):
                found = False
                for f in os.listdir(xqtl):
                    if fnmatch(f.replace('X','23').replace('Y','24'), f'*chr{chrom}.besd'):
                        xqtl_list.append(f'{xqtl}/'+f.replace('.besd',''))
                        found = True
                        break
                if found: continue
                # if there is no match
                xqtl_list.append(None)
        
        # parse input PLINK binaries
        if os.path.isfile(f'{bfile}.bed'):
            bfile_list = [bfile] * 24
        else:
            bfile_list = []
            for chrom in range(1,25):
                found = False
                for f in os.listdir(bfile):
                    if fnmatch(f.replace('X','23').replace('Y','24'), f'*chr{chrom}.bed'):
                        bfile_list.append(f'{bfile}/'+f.replace('.bed',''))
                        found = True
                        break
                if found: continue
                # if there is no match
                bfile_list.append(None)
        
        if os.path.isfile(f'{xqtl}.besd') and os.path.isfile(f'{bfile}.bed'):
            if not os.path.isfile(f'{out}/{prefix}.smr') or force:
                submitter.add(
                f'{smr} --bfile {bfile} --gwas-summary {tmpdir}/{prefix}.txt '+
                f'--beqtl-summary {xqtl} --out {out}/{prefix}.{qtl}'
                )
        
        else:
            # print(f'Processing: {prefix} \n\tConducting SMR by chromosome, following files have been found:')
            # for x, b in zip(xqtl_list, bfile_list):
            #     print('\t\t'.join([str(x), str(b)]))
            
            if not os.path.isdir(f'{out}/{prefix}.{qtl}'): os.mkdir(f'{out}/{prefix}.{qtl}')
            
            for x, b, chrom in zip(xqtl_list, bfile_list, range(1,25)):
                if x == None or b == None: continue
                if not os.path.isfile(f'{out}/{prefix}.{qtl}/chr{chrom}.smr') or force:
                    submitter.add(
                    f'{smr} --bfile {b} --gwas-summary {tmpdir}/{prefix}.txt '+
                    f'--beqtl-summary {x} --out {out}/{prefix}.{qtl}/chr{chrom}'
                    )
    
    qtl_list = []
    for y in os.listdir(args.qtl):
        if fnmatch(y, '*.besd'): qtl_list.append(f'{args.qtl}/{y}')
        if os.path.isdir(f'{args.qtl}/{y}'):
            if any([fnmatch(z, '*.besd') for z in os.listdir(f'{args.qtl}/{y}')]):
                qtl_list.append(f'{args.qtl}/{y}')
    
    print('Following QTL have been found:')
    for x in qtl_list: print(x)
    
    for x in args.pheno:
        if not os.path.isdir(f'{args.out}/{x}'): os.system(f'mkdir -p {args.out}/{x}')
        
        # scans directory for fastGWA files
        flist = []
        for y in os.listdir(f'{args._in}/{x}'):
            if fnmatch(y,'*.fastGWA') and (not fnmatch(y,'*_X.fastGWA')) and (not fnmatch(y, '*_all_chrs*')):
                flist.append(y)

        for y in flist:
            # conduct SMR for each QTL file
            for qtl in qtl_list:
                annot_smr(
                    gwa = f'{args._in}/{x}/{y}',
                    xqtl = qtl,
                    bfile = args.bfile,
                    out = f'{args.out}/{x}',
                    smr = args.smr,
                    force = args.force
                    )
    
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
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/smr') # intentionally absolute
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