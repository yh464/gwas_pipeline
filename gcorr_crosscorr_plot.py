#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2023-07-18
Version 2: 2024-11-14

A simplified script to plot genetic correlation between groups of phenotypes

Requires following inputs: 
    LDSC rg logs
Changelog:
    Changed the file name structure (now scans sumstats directory to find relevant rg log files)
    Changed the plotting style to scatterplot-style heatmaps
    Applied a wide- and long-format tabular output
'''

def main(args):
    import os
    from fnmatch import fnmatch
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import scipy.stats as sts
    
    # Red-blue colour map
    cdict = dict(red = ((0,0,0),(1/2,1,1),(1,.8,.8)),
                 green = ((0,0,0),(1/2,1,1),(1,0,0)),
                 blue = ((0,.8,.8),(1/2,1,1),(1,0,0)))
    cmap_name = 'redblue'
    cmap = mpl.colors.LinearSegmentedColormap(cmap_name,cdict,1024)
    try:
      mpl.colormaps.register(cmap)
    except:
      pass
    sns.set_theme(style = 'whitegrid')
    
    # scans directories to include sumstats 
    os.chdir(args.sumstats)
    prefix_1 = []; pheno_1 = []; count_1 = []
    prefix_2 = []; pheno_2 = []; count_2 = []
    for p in args.p1:
        c = 0
        for x in sorted(os.listdir(p)):
            if fnmatch(x,'*.sumstats'):
                prefix_1.append(x.replace('.sumstats','')); pheno_1.append(p)
                c += 1
        count_1.append(c)
    for p in args.p2:
        c = 0
        for x in sorted(os.listdir(p)):
            if fnmatch(x,'*.sumstats'):
                prefix_2.append(x.replace('.sumstats','')); pheno_2.append(p)
                c += 1
        count_2.append(c)

    # parse LDSC log files
    summary = []
    os.chdir(args._in)
    for p1,x1 in zip(pheno_1, prefix_1): # usually imaging phenotypes
        for p2,x2 in zip(pheno_2, prefix_2):
            fname = f'{p1}_{x1}.{p2}_{x2}.rg.log'
            if not os.path.isfile(fname):
                fname = f'{p2}_{x2}.{p1}_{x1}.rg.log'
            
            # The following section reads the RG and SE statistics from the log file
            if os.path.isfile(fname):
                tmp = open(fname)
                tmp_stats = tmp.read().splitlines()
                tmp_stats = tmp_stats[-4].split()
                while tmp_stats.count('') > 0:
                  tmp_stats.remove('')
                try: 
                    rg = float(tmp_stats[2])
                    if rg > 1: rg = 1
                    if rg < -1: rg = -1
                except: 
                  rg = np.nan
                  print(f'{fname} shows NA correlation!')
                try: se = max((float(tmp_stats[3]),10**-20))
                except: se = np.nan
            
                p = 1-sts.chi2.cdf((rg/se)**2, df = 1) # p value
                
                summary.append(pd.DataFrame(dict(group1 = p1, pheno1 = x1,
                  group2 = p2, pheno2 = x2, rg = rg, se = se, p = [p])))
    
    # statistics for the summary table
    summary = pd.concat(summary) # creates a long format table
    summary.insert(loc = len(summary.columns), column = 'q', value = np.nan)
    for p1,x1 in zip(pheno_1, prefix_1): # FDR correction for each IDP, which are non-independent
        summary.loc[(summary.pheno1==x1) & (summary.group1==p1) & ~np.isnan(summary.p),'q'] \
            = sts.false_discovery_control(
            summary.loc[(summary.pheno1==x1)&(summary.group1==p1) & ~np.isnan(summary.p),'p'])
    
    summary['pheno1'] = summary.pheno1.str.replace('_0.01','').replace('_meta','').replace('_corr','')
    summary['pheno2'] = summary.pheno2.str.replace('_meta','').replace('_0.01','')
    
    # tabular output, wide and long
    fout = f'{args.out}/crosscorr_' + '_'.join(args.p1)+'.'+'_'.join(args.p2)
    rg_tbl = summary.pivot(index = ['group1','pheno1'], columns = ['group2','pheno2'], values = 'rg')
    rg_tbl.to_csv(f'{fout}.wide.txt', index_label = False, sep = '\t',
                  header = True, index = True)
    summary.to_csv(f'{fout}.txt', index = False, sep = '\t')
    
    # plot figure
    from _plots import corr_heatmap
    fig = corr_heatmap(summary)
    fig.savefig(f'{fout}.pdf', bbox_inches = 'tight')
    
    ## relative widths/heights of plots are given in the 'count' section
    # fig, ax = plt.subplots(len(args.p1), len(args.p2), 
    #                      figsize = (sum(count_2)/3, sum(count_1)/3),
    #                      width_ratios = count_2, height_ratios = count_1)
    # if len(args.p1) == 1 and len(args.p2) == 1:
    #     ax = np.array([ax])
    # ax = ax.reshape((len(args.p1),len(args.p2)))
    # summary['pt_size'] = (summary.q < 0.05).astype(float) + \
    #     (summary.p < 0.05).astype(float) + 1 # FDR < 0.05 -> size = 3
    
    # for i in range(len(args.p1)):
    #     for j in range(len(args.p2)):
    #         p1 = args.p1[i]; p2 = args.p2[j]
    #         tmp = summary.loc[(summary.group1 == p1) & (summary.group2 == p2),:]
            
    #         sns.scatterplot(
    #             tmp,
    #             x = 'pheno2', y = 'pheno1',
    #             hue = 'rg', palette = 'redblue', hue_norm = (-1,1),
    #             size = 'pt_size', sizes = (25, 250),
    #             edgecolor = '.7',
    #             legend = False, ax = ax[i,j]
    #             )
            
    #         # remove x and y labels if aligned
    #         if i != len(args.p1) -1: 
    #             ax[i,j].set_xlabel(''); 
    #             ax[i,j].set_xticklabels([''] * len(ax[i,j].get_xticklabels()))
    #         else:
    #             ax[i,j].set_xlabel(args.p2[j])
    #             for label in ax[i,j].get_xticklabels():
    #               label.set_rotation(90)
    #         if j > 0:
    #             ax[i,j].set_ylabel(''); 
    #             ax[i,j].set_yticklabels([''] * len(ax[i,j].get_yticklabels()))
    #         else:
    #             ax[i,j].set_ylabel(args.p1[i])
            
    #         ax[i,j].set_xlim(-0.5, count_2[j]-0.5)
    #         ax[i,j].set_ylim(-0.5, count_1[i]-0.5)
            
    #         # despine
    #         for _, spine in ax[i,j].spines.items():
    #             spine.set_visible(False)
    # # colour bar
    # norm = mpl.colors.Normalize(vmin=-1, vmax=1)
    # plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap='redblue'), cax = fig.add_axes((0.92, 0.25, 0.02, 0.50)))
    # plt.savefig(f'{fout}.pdf', bbox_inches = 'tight')
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'This programme parses genetic cross-correlation between two groups of phenotypes')
    parser.add_argument('-p1', help = 'First group of phenotypes to correlate, usually IDP', nargs = '*')
    parser.add_argument('-p2', help = 'Second group of phenotypes to correlate, usually disorders', nargs = '*',
      default = ['global_structural','disorders','gradients'])
    parser.add_argument('-i','--in', dest = '_in', help = 'Directory for rg logs',
      default = '../gene_corr/gcorr/')
    parser.add_argument('--sumstats', help = 'sumstats directory to be scanned for file names',
      default = '../gene_corr/ldsc_sumstats/')
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory')
    # always forces output
    args = parser.parse_args()
    
    import os
    args._in = os.path.realpath(args._in)
    args.sumstats = os.path.realpath(args.sumstats)
    if type(args.out) == type(None): args.out = os.path.realpath(f'{args._in}/../')
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheno_%maf.%pheno.rg.log', __file__)
    proj.add_output(args.out+'/crosscorr_.*.pdf', __file__)
    try: main(args)
    except: cmdhistory.errlog()