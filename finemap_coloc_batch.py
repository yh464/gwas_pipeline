#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2024-11-25

A script to batch run hyprcoloc for groups of phenotypes

Preceding workflow:
    gwa_batch.py
Required input:
    GWAS summary stats for a group of traits
    clump output (sentinel variants)
'''

def main(args):
    import os
    import pandas as pd
    from fnmatch import fnmatch
    import numpy as np
    
    # array submitter
    from _utils import array_submitter
    submitter = array_submitter.array_submitter(
        name = 'coloc_'+'_'.join(args.pheno),
        env = 'gentoolsr',
        n_cpu = 1,
        timeout = 360,
        debug = True
        )
    
    from time import perf_counter as t
    tic = t()
    force = '-f' if args.force else ''
    
    # identify blocks of fine-mapping segments
    blocks = []
    # read clump outputs for each phenotype group
    for x in args.pheno:
        print(f'Identifying blocks for {x}')
        try:
            overlaps = pd.read_table(f'{args.clump}/{x}_{args.p:.0e}_overlaps.txt').set_index('label')
        except: # identify overlaps with lowest p-value above the p threshold
            flist = [] 
            for y in os.listdir(args.clump):
                if fnmatch(y,f'{x}_?e-??_overlaps.txt'): flist.append(y)
            plist = [float(z.split('_')[-2]) for z in flist]
            overlaps = pd.read_table(f'{args.clump}/{x}_{min(plist):.0e}_overlaps.txt').set_index('label')
            
        snps = overlaps.columns.tolist()
        
        # read a fastGWA file to determine the CHR and POS of the SNPs
        for y in os.listdir(f'{args._in}/{x}'):
            if fnmatch(y, '*.fastGWA') and not fnmatch(y, '*all_chrs*') and not fnmatch(y, '*_X.fastGWA'):
                tmp = pd.read_table(f'{args._in}/{x}/{y}', usecols = ['CHR','SNP','POS'], index_col = 'SNP')
                break

        # select 1MB chunks for fine-mapping
        tmp = tmp[['CHR','POS']]
        tmp.columns = ['chr','pos']
        tmp = tmp.loc[tmp.index.isin(snps), :]
        tmp['start'] = tmp['pos']- 5e+5
        tmp['stop'] = tmp['pos'] + 5e+5
        blocks.append(tmp[['chr','start','stop']])
        del tmp
    
    # identify overlaps between blocks from different phenotype groups
    blocks = pd.concat(blocks).sort_values(by = ['chr','start']).reset_index(drop=True)
    for x in range(1, blocks.shape[0]):
        if blocks.loc[x, 'start'] < blocks.loc[x-1, 'stop'] and \
            blocks.loc[x,'chr'] == blocks.loc[x-1,'chr']:
            blocks.loc[x,'start'] = blocks.loc[x-1, 'start']
            blocks.loc[x-1,['start','stop']] = np.nan
    blocks = blocks.dropna().reset_index(drop=True)
    toc = t()-tic
    print(f'Identified {blocks.shape[0]} blocks for multivariate fine-mapping, time = {toc:.3f}')
    
    # scans directory for fastGWA files
    flist = []
    for x in args.pheno:
        for y in os.listdir(f'{args._in}/{x}'):
            if fnmatch(y, '*.fastGWA') and not fnmatch(y, '*all_chrs*') and not fnmatch(y, '*_X.fastGWA'):
                flist.append(f'{args._in}/{x}/{y}')
    print(f'Found {len(flist)} GWAS summary statistics files.')
    
    # specify output
    outdir = f'{args.out}/'+'_'.join(sorted(args.pheno))
    if not os.path.isdir(outdir): os.system(f'mkdir -p {outdir}')
    for x in range(1, blocks.shape[0]):
        c = blocks.loc[x, 'chr']
        start = blocks.loc[x,'start']
        stop = blocks.loc[x,'stop']
        out = f'{outdir}/chr{c:.0f}_{start:.0f}_{stop:.0f}_coloc.txt'
        if not os.path.isfile(out) or args.force:
            cmd = 'Rscript finemap_hyprcoloc.r -i '+':'.join(flist) + \
                f' --chr {c:.0f} --start {start:.0f} --stop {stop:.0f} -o {out} {force}'
            submitter.add(cmd)
    submitter.submit()
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
      description = 'This programme batch runs the fine-map pipeline')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
    parser.add_argument('-i','--in', dest = '_in', help = 'Directory containing all summary stats',
      default = '../gwa/')
    parser.add_argument('-c','--clump', dest = 'clump', help = 'Directory containing all clump outputs',
      default = '../clump/')
    parser.add_argument('-o', '--out', dest = 'out', help = 'output directory',
      default = '../coloc/')
    parser.add_argument('-p', dest = 'p', help = 'p-value', default = 3.1076e-11, type = float)
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
      default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','clump']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
    from _utils import path, cmdhistory, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_var('%pval',r'[0-9.+-e]+', 'minor allele freq') # only allows digits and decimals
    proj.add_input(args._in+'/%pheng/%pheno_%maf.fastGWA', __file__)
    proj.add_input(args.clump+'/%pheng_%pval_overlaps.txt',__file__)
    proj.add_output(args.out+'/%pheng/*',__file__)
    try: main(args)
    except: cmdhistory.errlog()