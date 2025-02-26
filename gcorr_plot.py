#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2024-11-18

Plots genetic correlation within groups of phenotypes, and heritability

Requires following inputs: 
    LDSC rg logs
Upstream input:
    heri_batch.py
    gcorr_batch.py
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
    from _utils.path import normaliser
    
    normer = normaliser()
    
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
    prefix_list = []; pheno_list = []; counts = []
    for p in args.pheno:
        c = 0
        for x in sorted(os.listdir(p)):
            if not fnmatch(x,'*.sumstats'): continue
            if any([x.find(ex) > -1 for ex in args.exclude]): continue
            prefix_list.append(x.replace('.sumstats','')); 
            pheno_list.append(p)
            c += 1
        counts.append(c)
    
    # parse LDSC log files
    summary = []
    os.chdir(args._in)
    for i in range(len(pheno_list)):
        p1 = pheno_list[i]; x1 = prefix_list[i]
        
        # h2 logs
        fname = f'{args.sumstats}/{p1}/{x1}.h2.log'
        if os.path.isfile(fname):
            tmp = open(fname).read().splitlines()[-7].replace('(','').replace(')','').split()
            while tmp.count('') > 0:
                tmp.remove('')
            try: h2 = float(tmp[-2])
            except: h2 = np.nan
            try: se = float(tmp[-1])
            except: se = np.nan
            
            p = 1-sts.chi2.cdf((h2/se)**2, df = 1)
            summary.append(pd.DataFrame(dict(group1 = p1, pheno1 = x1,
              group2 = p1, pheno2 = x1, rg = h2, se = se, p = [p]))) # we hack the table a bit, pretend as if rg is h2
            
        # rg logs
        for j in range(i):
            p2 = pheno_list[j]; x2 = prefix_list[j]
            if p1 < p2: fname = f'{p1}.{p2}/{p1}_{x1}.{p2}_{x2}.rg.log'
            else: fname = f'{p2}.{p1}/{p2}_{x2}.{p1}_{x1}.rg.log'
            
            # check for irregularly named files:
            for fname1 in [f'{p1}.{p2}/{p2}_{x2}.{p1}_{x1}.rg.log',
                           f'{p1}.{p2}/{p1}_{x1}.{p2}_{x2}.rg.log',
                           f'{p2}.{p1}/{p2}_{x2}.{p1}_{x1}.rg.log',
                           f'{p2}.{p1}/{p1}_{x1}.{p2}_{x2}.rg.log']:
                if os.path.isfile(fname1) and not os.path.isfile(fname):
                    os.rename(fname1, fname)
                if os.path.isfile(fname1) and os.path.isfile(fname) and fname1 != fname:
                    os.remove(fname1)
            
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
    summary['q'] = np.nan
    summary.loc[~np.isnan(summary.p),'q'] = sts.false_discovery_control(summary.loc[~np.isnan(summary.p),'p'])
    summary['pheno1'] = summary.pheno1.str.replace('_0.01','').str.replace('_meta','').str.replace('_corr','')
    summary['pheno2'] = summary.pheno2.str.replace('_0.01','').str.replace('_meta','').str.replace('_corr','')
    
    # tabular output, wide and long
    fout = f'{args.out}/corr_' + '_'.join(args.pheno)
    rg_tbl = summary.pivot(index = ['group1','pheno1'], columns = ['group2','pheno2'], values = 'rg')
    
    # flip the table around for plotting
    summary1 = summary[['group2','pheno2','group1','pheno1','rg','se','p','q']]
    summary1.columns = summary.columns
    summary = pd.concat([summary, summary1], axis = 0).drop_duplicates()
    
    # plot figure
    ## relative widths/heights of plots are given in the 'count' section
    fig, ax = plt.subplots(len(args.pheno), len(args.pheno), 
                         figsize = (sum(counts)/3, sum(counts)/3),
                         width_ratios = counts, height_ratios = counts[::-1])
    try: ax = ax.reshape((len(args.pheno),len(args.pheno)))
    except: ax = np.array([ax]).reshape((len(args.pheno),len(args.pheno)))
    summary['pt_size'] = (summary.q < 0.05).astype(float) + \
        (summary.p < 0.05).astype(float) + 1 # FDR < 0.05 -> size = 3
    
    for i in range(len(args.pheno)):
        for j in range(len(args.pheno)):
            p1 = args.pheno[-1-i]; p2 = args.pheno[j]
            tmp = summary.loc[(summary.group1 == p1) & (summary.group2 == p2),:]
            
            sns.scatterplot(
                tmp,
                x = 'pheno2', y = 'pheno1',
                hue = 'rg', palette = 'redblue', hue_norm = (-1,1),
                size = 'pt_size', sizes = (25, 250), size_norm = (1,3),
                edgecolor = '.7',
                legend = False, ax = ax[i,j]
                )
            
            # remove x and y labels if aligned
            if i != len(args.pheno) -1: 
                ax[i,j].set_xlabel(''); 
                ax[i,j].set_xticklabels([''] * len(ax[i,j].get_xticklabels()))
            else:
                ax[i,j].set_xlabel(args.pheno[j])
                for label in ax[i,j].get_xticklabels():
                  label.set_rotation(90)
            if j > 0:
                ax[i,j].set_ylabel(''); 
                ax[i,j].set_yticklabels([''] * len(ax[i,j].get_yticklabels()))
            else:
                ax[i,j].set_ylabel(args.pheno[-1-i])
            
            ax[i,j].set_xlim(-0.5, counts[j]-0.5)
            ax[i,j].set_ylim(-0.5, counts[-1-i]-0.5)
            
            # despine
            for _, spine in ax[i,j].spines.items():
                spine.set_visible(False)
            
            # diagonal lines
            if args.pheno[-1-i] == args.pheno[j]:
                ax[i,j].axline((0,0), slope = 1, color = 'k', zorder = 0)
            
    ax[-1,0].annotate('Heritability', (0,0),xytext=(-4,-3.5), rotation = 45)

    # colour bar
    norm = mpl.colors.Normalize(vmin=-1, vmax=1)
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap='redblue'), cax = fig.add_axes((0.92, 0.25, 0.02, 0.50)))
    
    # output
    plt.savefig(f'{fout}.pdf', bbox_inches = 'tight')
    normer.normalise(summary).to_csv(f'{fout}.txt', index = False, sep = '\t')
    normer.normalise(rg_tbl).to_csv(f'{fout}.wide.txt', index_label = False, sep = '\t',
                  header = True, index = True)
    
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'This programme parses genetic cross-correlation between two groups of phenotypes')
    parser.add_argument('pheno', help = 'Phenotypes to correlate', nargs = '*')
    parser.add_argument('-i','--in', dest = '_in', help = 'Directory for rg logs',
      default = '../gcorr/rglog/')
    parser.add_argument('--sumstats', help = 'sumstats directory to be scanned for file names',
      default = '../gcorr/ldsc_sumstats/')
    parser.add_argument('--exclude', help = 'phenotypes to exclude', nargs = '*', default = [])
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