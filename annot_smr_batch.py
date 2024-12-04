#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2024-11-08

Conducts Summary-data Mendelian Randomisation to prioritise genes using
expression, isoform, splicing, methylation and chromosome accessibility QTL.

Requires following inputs: 
    fine-map extracted list of SNPs, (finemap_batch.py, finemap_parse.py)
    GWAS summary statistics,
    BESD format xQTL datasets.
'''

def format_gwa(gwa, snp_list, tmpgwa):
    # input file name, usually fastGWA format
    import pandas as pd
    df = pd.read_table(gwa, sep = '\s+')
    df = pd.merge(df, snp_list)
    if not 'N' in df.columns:
        df['N'] = 'NA' # will not be used in SMR, just for the correct format 
    df = df[['SNP','A1','A2','AF1','BETA','SE','P','N']]
    df.to_csv(tmpgwa, sep = '\t', header = True, index = False)
    
def main(args):
    import os
    from fnmatch import fnmatch
    import pandas as pd
    
    # temp directory
    tmpdir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/temp/smr_cache'
    if not os.path.isdir(tmpdir): os.system(f'mkdir -p {tmpdir}')
    
    # array submitter
    from _utils import array_submitter
    submitter = array_submitter.array_submitter(
        name = f'genpr_smr_{args.pheno}',
        timeout = 90,
        debug = True
    )
    
    qtl_list = []
    for y in os.listdir(args.qtl):
        if fnmatch(y, '*.besd'): qtl_list.append(y)
    
    for x in args.pheno:
        if not os.path.isdir(f'{args.out}/{x}'): os.system(f'mkdir -p {args.out}/{x}')
        
        # list of sig. variants
        sig_snp = args.finemap.replace('%pheno',x)
        sig_snp = pd.read_table(sig_snp, header = None)
        sig_snp.columns = ['SNP']
        
        # scans directory for fastGWA files
        flist = []
        for y in os.listdir(f'{args._in}/{x}'):
            if fnmatch(y,'*.fastGWA') and (not fnmatch(y,'*_X.fastGWA')) and (not fnmatch(y, '*_all_chrs*')):
                flist.append(y)

        for y in flist:
            # format GWAS
            tmpgwa = f'{tmpdir}/{x}/{y}'
            if args.force or not os.path.isfile(tmpgwa):
                format_gwa(y, sig_snp, tmpgwa)
            
            # conduct SMR for each QTL file
            for qtl in qtl_list:
                out_prefix = f'{args.out}/{x}/{y}_{qtl}'.replace('.fastGWA','').replace('.besd','')
                if os.path.isfile(f'{out_prefix}.smr') and not args.force: continue
                submitter.add(f'{args.smr} --bfile {args.bfile} --gwas-summary {args._in}/{x}/{y} '+
                              f'--beqtl-summary {args.qtl}/{qtl} --out {out_prefix} --thread-num 4')
    
    submitter.submit()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
      description = 'This programme batch runs summary data randomisation')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
    parser.add_argument('-i','--in', dest = '_in', help = 'Directory containing all GWA summary statistics',
      default = '../gwa/')
    parser.add_argument('-q','--qtl', dest = 'qtl', help = 'Directory containing all xQTL files',
      default = '../params/xqtl')
    parser.add_argument('-s','--smr', dest = 'smr', help = 'Location of SMR binary',
      default = '/rds/project/rds-Nl99R8pHODQ/UKB/Imaging_genetics/yh464/toolbox/smr') # intentionally absolute
    parser.add_argument('-f', '--finemap', dest = 'finemap', 
      help = 'file name of fine-mapped SNP lists, %pheno for phenotype wildcard',
      default = '../finemap/%pheno_sig_variants_list.txt')
    parser.add_argument('-b', '--bfile', dest = 'bfile', help = 'bed binary to use in magma',
      default = '/rds/project/rds-Nl99R8pHODQ/UKB/Imaging_genetics/yh464/bed/autosomes') # intentionally absolute
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
      default = '../genpr/smr')
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
      default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','finemap','smr','bfile']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheng/%pheno_%maf.fastGWA', __file__)
    proj.add_output(args.out+'/%pheng/%pheno_%maf.finemap.summary',__file__)
    try: main(args)
    except: cmdhistory.errlog()