#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-01-07

A generalised plotting tool to plot scatterplot-style heatmaps for correlations
'''

def corr_heatmap(summary, absmax = None):
    '''
    Required input format: long format pd.DataFrame
        1st column: group label (x axis, as xlabel)
        2nd column: individual label (x axis, as xticklabels)
        3rd column: group label (y axis, as ylabel)
        4th column: individual label (y axis, as yticklabels)
        5th column: correlation value
        and may contain a 'p' and 'q' column somewhere in the data frame
    '''
    import seaborn as sns
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import numpy as np
    from scipy.stats import false_discovery_control as fdr
    
    # Red-blue colour map
    cdict = dict(red = ((0,0,0),(1/2,1,1),(1,.8,.8)),
                 green = ((0,0,0),(1/2,1,1),(1,0,0)),
                 blue = ((0,.8,.8),(1/2,1,1),(1,0,0)))
    cmap_name = 'redblue'
    cmap = mpl.colors.LinearSegmentedColormap(cmap_name,cdict,1024)
    try: mpl.colormaps.register(cmap)
    except: pass
    
    # style sheet
    sns.set_theme(style = 'whitegrid')
    
    # determine figure size and aspect ratios
    group1 = summary.iloc[:, 0].unique()
    group2 = summary.iloc[:, 2].unique()
    
    count1 = [summary.loc[summary.iloc[:,0] == x, summary.columns[1]].unique().size for x in group1]
    count2 = [summary.loc[summary.iloc[:,2] == x, summary.columns[3]].unique().size for x in group2]
    
    fig, ax = plt.subplots(len(group1), len(group2),
                           figsize = (sum(count2)/3, sum(count1)/3),
                           width_ratios = count2, height_ratios = count1)
    
    if len(group1) == 1 and len(group2) == 1: ax = np.array([ax])
    ax = ax.reshape((len(group1), len(group2)))
    
    # estimate significance:
    if 'fdr' in summary.columns: summary['q'] = summary['fdr']
    if 'FDR' in summary.columns: summary['q'] = summary['fdr']
    if 'P' in summary.columns: summary['p'] = summary['P']
    
    if not 'q' in summary.columns and 'p' in summary.columns:
        summary['q'] = np.nan
        for g in group1:
            for trait in summary.loc[summary.iloc[:,0] == g, summary.columns[1]]:
                summary.loc[(summary.iloc[:,0] == g) & (summary.iloc[:,1] == trait) & ~np.isnan(summary['p']),'q'] = \
                fdr(summary.loc[(summary.iloc[:,0] == g) & (summary.iloc[:,1] == trait) & ~np.isnan(summary['p']),'p'])
    
    if 'p' in summary.columns:
        summary['Significance'] = 'NS'
        summary.loc[summary['p'] < 0.05, 'Significance'] = 'nominal'
        summary.loc[summary['q'] < 0.05, 'Significance'] = 'FDR-sig'
        legend_param = 'brief'
    else:
        summary['Significance'] = 'NA'
        legend_param = False
    # point size will depend on significance
    sizes = {'NS': 25, 'nominal': 125, 'FDR-sig': 250, 'NA': 250}
    
    # colour bar range (may manually set to 1)
    if type(absmax) == type(None): absmax = np.abs(summary.iloc[:,4]).max()
    
    for i in range(len(group1)):
        for j in range(len(group2)):
            g1 = group1[i]; g2 = group2[j]
            tmp = summary.loc[(summary.iloc[:,0] == g1) & (summary.iloc[:,2] == g2),:]
            sns.scatterplot(
                tmp,
                x = summary.columns[3], y = summary.columns[1],
                hue = summary.columns[4], hue_norm = (-absmax, absmax),
                palette = 'redblue',
                size = 'Significance', sizes = sizes,
                edgecolor = '.7',
                # legend = legend_param if i == len(group1) - 1 and j == 0 else False,
                legend = False,
                ax = ax[i,j]
                )
            
            # remove x labels except for last row
            if i != len(group1) - 1:
                ax[i,j].set_xlabel('')
                ax[i,j].set_xticklabels([''] * len(ax[i,j].get_xticklabels()))
            else:
                ax[i,j].set_xlabel(group2[j])
                for label in ax[i,j].get_xticklabels():
                    label.set_rotation(90)
            
            # remove y labels except for first column
            if j > 0:
                ax[i,j].set_ylabel(''); 
                ax[i,j].set_yticklabels([''] * len(ax[i,j].get_yticklabels()))
            else:
                ax[i,j].set_ylabel(group1[i])
            
            # set x and y limits
            ax[i,j].set_xlim(-0.5, count2[j]-0.5)
            ax[i,j].set_ylim(-0.5, count1[i]-0.5)
            
            # despine
            for _, spine in ax[i,j].spines.items():
                spine.set_visible(False)
    
    # colour bar
    norm = mpl.colors.Normalize(vmin = -absmax, vmax = absmax)
    cax = fig.add_axes((0.95, 0.4, 0.03, 0.35))
    cax.set_title(summary.columns[4])
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap='redblue'), cax = cax)
    
    # move legend location
    if legend_param:
        # leg = ax[len(group1)-1, 0].get_legend()
        # handles = leg.legend_handles
        # remove colour bar related items
        # handles = [handles[i] for i in [-3, -1, -2]]
        # ax[len(group1)-1,0].legend().remove()
        
        handles = [
            mpl.lines.Line2D([],[], marker = 'o', linewidth = 0, markersize = 25**0.5, # Seaborn point size correspond to the square of diameter 
                markeredgecolor = '.7', markerfacecolor = '.7', label = 'not significant'),
            mpl.lines.Line2D([],[], marker = 'o', linewidth = 0, markersize = 125**0.5, 
                markeredgecolor = '.7', markerfacecolor = '.7', label = 'nominal'),
            mpl.lines.Line2D([],[], marker = 'o', linewidth = 0, markersize = 250**0.5, 
                markeredgecolor = '.7', markerfacecolor = '.7', label = 'FDR-corrected')
            ]
        fig.legend(handles=handles, title = 'Significance',
                   frameon = False, loc = 'upper left', bbox_to_anchor=(0.91, 0.4))
        
    # ax[len(group1)-1, 0].legend(bbox_to_anchor=(0.2, 1.05))
    
    return fig