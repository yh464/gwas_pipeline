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
    from hashlib import sha256
    from time import perf_counter as t
    tic = t()
    force = '-f' if args.force else ''
    tmpdir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/temp/coloc'
    if not os.path.isdir(tmpdir): os.system(f'mkdir -p {tmpdir}')

    # scans directory for fastGWA files
    from _utils.path import find_gwas
    from _plugins.logparser import crosscorr_parse
    gwa = []
    if len(args.pheno) > 1 and len(args.filter) > 0:
        if any([not x in args.pheno for x in args.filter]):
            raise ValueError('Filter must be based on traits included in the coloc analysis')
        gwa1 = find_gwas(args.filter, dirname = args._in, clump = True)
        other_traits = [x for x in args.pheno if x not in args.filter]
        gwa2 = find_gwas(other_traits, dirname = args._in, clump = True)
        rg = crosscorr_parse(gwa1, gwa2, logdir = args.rg)
        if args.rgp == None: rg = rg.loc[rg.q < 0.05, ['group1','pheno1','group2','pheno2']]
        else: rg = rg.loc[rg.p < args.rgp, ['group1','pheno1','group2','pheno2']]
        
        # include phenotypes from the --filter group if correlated with anything else
        # include phenotypes from other groups if correlated with anything in --filter
        rg1 = rg[['group1','pheno1']]; rg2 = rg[['group2','pheno2']]
        rg2.columns = ['group1','pheno1']
        rg = pd.concat([rg1,rg2]).drop_duplicates().reset_index(drop = True)
        for idx in rg.index:
            x,y = rg.loc[idx,:]
            gwa.append((x,y))
    else:
        gwa = find_gwas(*args.pheno, dirname = args._in, long = True)
    print(f'Found {len(gwa)} GWAS summary statistics files.')
    print(gwa)
    
    # identify blocks of fine-mapping segments
    from _plugins.logparser import parse_clump
    _, loci = parse_clump(gwa, clump_dir = args.clump, pval = args.pval)
    loci = loci.loc[loci.P < args.pval, ['CHR', 'START', 'STOP']]
    loci['START'] -= 5e5; loci['STOP'] += 5e5
    loci = loci.dropna().reset_index(drop=True)
    toc = t()-tic
    print(f'Identified {loci.shape[0]} blocks for multivariate fine-mapping, time = {toc:.3f}')
    
    # check cache
    cache_files = [f'{args.out}/loci/{g}/{p}_chr{loci.CHR.iloc[i]}_{loci.START.iloc[i]:.0f}_{loci.STOP.iloc[i]:.0f}.txt'\
        for g,p in gwa for i in range(loci.shape[0])]
    if not all([os.path.isfile(x) for x in cache_files]):
        import finemap_coloc_extract
        dependency = finemap_coloc_extract.main(
            pheno = [f'{g}/{p}' for g,p in gwa],
            _in = args._in, clump = args.clump, pval = args.pval, out = args.out,
            debug = args.debug
        )
    else: dependency = []

    # array submitter
    from _utils.slurm import array_submitter
    submitter = array_submitter(
        name = 'coloc_'+'_'.join(args.pheno),
        env = 'gentoolsr', n_cpu = 1,
        timeout = len(gwa)+10, dependency = dependency
        )
    
    # correlation matrix is needed for flashfm
    if args.flashfm:
        rg = crosscorr_parse(gwa)
        rg['p1'] = rg.group1 + '_' + rg.pheno1; rg['p2'] = rg.group2 + '_' + rg.pheno2
        rg = rg.pivot_table(index = 'p1', columns = 'p2', values = 'rg')
        rg.columns.name = None; rg.index.name = None
        rg = rg.fillna(rg.T)
        for x in rg.columns: rg.loc[x,x] = 1
        rg.to_csv(f'{tmpdir}/{sha256(repr(gwa).encode()).hexdigest()[:10]}_rg.txt', sep = '\t')

    outdir = f'{args.out}/'+'_'.join([g for g,_ in find_gwas(args.pheno)])
    if not os.path.isdir(outdir): os.system(f'mkdir -p {outdir}')
    for x in range(1, loci.shape[0]):
        c = loci.loc[x, 'CHR']
        start = loci.loc[x,'START']
        stop = loci.loc[x,'STOP']
        if args.hyprcoloc:
            out = f'{outdir}/chr{c:.0f}_{start:.0f}_{stop:.0f}_hyprcoloc.txt'
            if not os.path.isfile(out) or args.force:
                cmd = ['Rscript','finemap_hyprcoloc.r'] + [f'{g}/{p}' for g,p in gwa] + \
                    [f'-i {args.out}/loci --chr {c:.0f} --start {start:.0f} --stop {stop:.0f} -o {out} {force}']
                submitter.add(' '.join(cmd))
        if args.flashfm:
            out = f'{outdir}/chr{c:.0f}_{start:.0f}_{stop:.0f}_flashfm'
            if not os.path.isfile(out+'.txt') or args.force:
                cmd = ['Rscript', 'finemap_flashfm.r'] + [f'{g}/{p}' for g,p in gwa] + \
                    ['-i', f'{args.out}/loci',f'--chr {c:.0f} --start {start:.0f} --stop {stop:.0f} -o {out}',
                    '--gcov',f'{tmpdir}/{sha256(repr(gwa).encode()).hexdigest()[:10]}_rg.txt', force]
                submitter.add(' '.join(cmd))
    submitter.submit()
    
if __name__ == '__main__':
    from _utils.slurm import slurm_parser
    parser = slurm_parser(
      description = 'This programme batch runs the fine-map pipeline')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
    parser.add_argument('-i','--in', dest = '_in', help = 'Directory containing all summary stats',
      default = '../gwa/')
    parser.add_argument('-c','--clump', dest = 'clump', help = 'Directory containing all clump outputs',
      default = '../clump/')
    parser.add_argument('-o', '--out', dest = 'out', help = 'output directory',
      default = '../coloc/')
    parser.add_argument('--filter', help = 'Filter for significant correlates of --filter phenotypes',
      nargs = '*', default = [])
    parser.add_argument('-r','--rg', help = 'Directory for rg logs, to filter traits',
      default = '../gcorr/rglog/')
    parser.add_argument('--rgp', help = 'p-value threshold for genetic correlation', default = None, type = float)
    parser.add_argument('--hyprcoloc', help = 'run hyprcoloc', action = 'store_true', default = False)
    parser.add_argument('--flashfm', help = 'run flashfm', action = 'store_true', default = False)
    parser.add_argument('--mvsusie', help = 'run mvsusie', action = 'store_true', default = False)
    parser.add_argument('-p', '--pval', help = 'p-value', default = 5e-8, type = float)
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
      default = False, action = 'store_true')
    args = parser.parse_args()

    # pre-process cmd line arguments
    import os
    for arg in ['_in','out','clump', 'rg']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    if args.rgp != None and 0 < args.rgp < 1:
        if len(args.filter) == 0: Warning('Ignoring rgp value as no filtering trait is specified')
    if args.rgp != None and (args.rgp > 1 or args.rgp <= 0): raise ValueError('p-value threshold must be 0 to 1')
    if not any([args.hyprcoloc, args.flashfm, args.mvsusie]):
        Warning('No algorithm specified, defaulting to hyprcoloc')
        args.hyprcoloc = True
    args.pheno.sort()

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