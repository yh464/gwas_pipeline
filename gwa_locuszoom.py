#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-01-21

Plots locuszoom plots for any number of input GWAS summary stats and one single SNP
'''

def main(args):
    import os
    # output directory
    if not os.path.isdir(os.path.dirname(args.out)): 
        os.system(f'mkdir -p {os.path.dirname(args.out)}')
    
    # check progress
    if os.path.isfile(args.out) and not args.force: return
    
    # extract chunks containin index SNP
    all_sumstats = []
    import pandas as pd
    import numpy as np
    for file in args.gwa:
        prefix = '.'.join(os.path.basename(file).split('.')[:-1])
        prefix = prefix.replace('_0.01','').replace('_meta','')
        tmp = pd.read_table(file, usecols = ['CHR','SNP','POS','P','A1','A2'], index_col = 'SNP')
        
        if args.chr == None or args.pos == None:
            chrom = tmp.loc[args.snp, 'CHR']; pos = tmp.loc[args.snp, 'POS']
        else:
            chrom = args.chr; pos = args.pos
        
        tmp = tmp.loc[(tmp.CHR == chrom) & (tmp.POS > pos - args.ld) & (tmp.POS < pos + args.ld), ['CHR','POS','A1','A2','P']]
        tmp['POS'] /= 1000000
        tmp['-log(P)'] = -np.log10(tmp.P)
        tmp.insert(loc = 0, column = 'Phenotype', value = prefix)
        all_sumstats.append(tmp)
    
    # plot locus zoom plot
    import seaborn as sns
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize = (7,3))
    
    # export sumstats
    all_sumstats = pd.concat(all_sumstats)
    all_sumstats.to_csv(args.out + '.txt', sep = '\t', index = False)
    
    sns.scatterplot(all_sumstats, x = 'POS', y = '-log(P)', hue = 'Phenotype', 
        palette = 'pastel', ax = ax, marker = '.', s = 24, linewidths = 0)
    ax.set_xlabel(f'Chromosome {chrom} (Mb)')
    ax.axhline(y = -np.log10(5e-8), c = 'k')
    
    # move legend location
    handles = ax.get_legend().legend_handles
    fig.legend(handles = handles, frameon = False, loc = 'center left', bbox_to_anchor = (0.91, 0.5))
    ax.legend().remove()
    
    fig.savefig(args.out, bbox_inches = 'tight')
    plt.close()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'This programme compiles LocusZoom plots for any number of input GWAS summary stats')
    parser.add_argument('gwa', help = 'GWAS summary stats', nargs = '*')
    parser.add_argument('-s','--snp', dest = 'snp', help = 'index SNP')
    parser.add_argument('-c','--chr', dest = 'chr', type = int, help = 'chromosome') # these take priority over SNP
    parser.add_argument('-p','--pos', dest = 'pos', type = int, help = 'genomic position') # these take priority over SNP
    parser.add_argument('-l', '--ld', help = 'flanking window', type = int, default = 1000000)
    parser.add_argument('-o','--out', dest = 'out', help = 'output file name')
    parser.add_argument('-f','--force',dest = 'force', help = 'force overwrite',
      default = False, action = 'store_true')
    args = parser.parse_args()

    # path manipulations
    import os
    args.out = os.path.realpath(args.out)
    if not args.out[-3:] in ['pdf','png','eps','svg']: args.out += '.pdf'
    args.gwa = [os.path.realpath(x) for x in args.gwa]
    args.gwa.sort()

    from _utils import logger
    logger.splash(args)
    main(args)