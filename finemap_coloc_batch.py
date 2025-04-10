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
        # debug = True
        )
    
    from time import perf_counter as t
    tic = t()
    force = '-f' if args.force else ''
    
    # scans directory for fastGWA files
    from _utils.path import find_gwas, find_clump
    flist = []
    gwa = []
    if len(args.pheno) > 1 and args.filter:
        gwa1 = find_gwas(args.pheno[0], dirname = args._in)
        gwa2 = find_gwas(*args.pheno[1:], dirname = args._in)
        from gcorr_plot import crosscorr_parse
        rg = crosscorr_parse(gwa1, gwa2, logdir = args.rg)
        rg = rg.loc[rg.q < 0.05, ['group1','pheno1','group2','pheno2']]
        
        # include phenotypes from group 1 if correlated with anything else
        # include phenotypes from other groups if correlated with anything in group 1
        rg1 = rg[['group1','pheno1']]; rg2 = rg[['group2','pheno2']]
        rg2.columns = ['group1','pheno1']
        rg = pd.concat([rg1,rg2]).drop_duplicates().reset_index(drop = True)
        for idx in rg.index:
            x,y = rg.loc[idx,:]
            flist.append(f'{args._in}/{x}/{y}.fastGWA')
            gwa.append((x,y))
    else:
        gwa = find_gwas(*args.pheno, dirname = args._in)
        for x, ys in gwa:
            for y in ys:
                flist.append(f'{args._in}/{x}/{y}.fastGWA')
                gwa.append((x,y))
    print(f'Found {len(flist)} GWAS summary statistics files.')
    submitter.config(timeout = len(flist) * 5)
    
    # identify blocks of fine-mapping segments
    blocks = []
    # read clump outputs for each phenotype
    for x,y in gwa:
        print(f'Identifying blocks for {x}/{y}')
        clump,_ = find_clump(f'{args.clump}/{x}',y,args.p)
        try: clump = pd.read_table(clump, sep = '\\s+', usecols = ['CHR','BP'])
        except: continue
        
        # select 1MB chunks for fine-mapping
        clump = clump[['CHR','BP']]
        clump.columns = ['chr','pos']
        clump['start'] = clump['pos']- 5e+5
        clump['stop'] = clump['pos'] + 5e+5
        blocks.append(clump[['chr','start','stop']])
    
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
    parser.add_argument('--filter', help = 'Filter for significantly correlated phenotypes',
      default = False, action = 'store_true')
    parser.add_argument('-r','--rg', help = 'Directory for rg logs, to filter traits',
      default = '../gcorr/rglog/')
    parser.add_argument('-p', dest = 'p', help = 'p-value', default = 5e-8, type = float)
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
      default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','clump', 'rg']:
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