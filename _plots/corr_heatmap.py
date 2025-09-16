#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-01-07
Version 2: 0225-04-03

A generalised plotting tool to plot scatterplot-style heatmaps for correlations
'''

def capitalise(series):
    out = []
    for x in range(len(series)):
        tmp = series.iloc[x]
        if type(tmp) != str: out.append(tmp)
        elif len(tmp) < 1: out.append(tmp)
        else:
            tmp = tmp[0].upper() + tmp[1:]
            out.append(tmp)
    return out

def corr_heatmap(summary, sort = True, absmax = None, autocor = False, annot = '', p_threshold: list[float] = []):
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
    import warnings
    from aes import redblue

    try: mpl.colormaps.register(redblue)
    except: pass
    
    # style sheet
    sns.set_theme(style = 'whitegrid')
    
    # normalise labels
    for col in range(4):
        summary.iloc[:,col] = summary.iloc[:,col].str.replace('_',' ').str.replace(
            ' l$',' L', regex = True).str.replace(' r$',' R', regex = True)
        summary.iloc[:,col] = capitalise(summary.iloc[:,col])
        
    # determine figure size and aspect ratios
    group1 = summary.iloc[:, 0].unique(); group1.sort()
    group2 = summary.iloc[:, 2].unique(); group2.sort()
    if len(group1) == len(group2):
        if all(group1 == group2): autocor = True # diagonal lines to be plotted if auto-correlating
    
    count1 = [summary.loc[summary.iloc[:,0] == x, summary.columns[1]].unique().size for x in group1]
    count2 = [summary.loc[summary.iloc[:,2] == x, summary.columns[3]].unique().size for x in group2]
    
    fig, ax = plt.subplots(len(group1), len(group2),
                           figsize = (sum(count2)/3, sum(count1)/3),
                           width_ratios = count2, height_ratios = count1[::-1], squeeze = False)
    ax = ax[::-1,:] # invert y axis
    
    # estimate significance:
    p_threshold = sorted(p_threshold, reverse = True)
    if 'fdr' in summary.columns: summary['q'] = summary['fdr']
    if 'FDR' in summary.columns: summary['q'] = summary['fdr']
    if 'P' in summary.columns: summary['p'] = summary['P']
    
    if not 'q' in summary.columns and 'p' in summary.columns:
        summary = summary.assign(q = np.nan)
        for g in group1:
            for trait in summary.loc[summary.iloc[:,0] == g, summary.columns[1]]:
                summary.loc[(summary.iloc[:,0] == g) & (summary.iloc[:,1] == trait) & ~np.isnan(summary['p']),'q'] = \
                fdr(summary.loc[(summary.iloc[:,0] == g) & (summary.iloc[:,1] == trait) & ~np.isnan(summary['p']),'p'])
    
    if 'p' in summary.columns and len(p_threshold) == 0:
        summary = summary.assign(Significance = 'not significant')
        summary.loc[summary['p'] < 0.05, 'Significance'] = 'nominal'
        summary.loc[summary['q'] < 0.05, 'Significance'] = 'FDR-corrected'
        sig_label = 'FDR-corrected'
        legend_param = 'brief'
        sizes = {'not significant': 25, 'nominal': 125, 'FDR-corrected': 250}
    elif 'p' in summary.columns and len(p_threshold) > 0:
        summary = summary.assign(Significance = 'not significant')
        for p_thr in p_threshold:
            summary.loc[summary['p'] < p_thr, 'Significance'] = f'p < {p_thr}'
        sig_label = f'p < {p_threshold[-1]}'
        legend_param = 'brief'
        sizes = dict(zip(['not significant'] + [f'p < {x}' for x in p_threshold], np.linspace(25, 250, len(p_threshold)+1)))
    else:
        summary = summary.assign(Significance = 'NA')
        sig_label = ''
        legend_param = False
        sizes = {'NA': 250}
    
    # colour bar range (may manually set to 1)
    if type(absmax) == type(None): absmax = np.abs(summary.iloc[:,4]).max()
    
    # order of columns and rows
    if not sort:
        map_order = dict(zip(summary.iloc[:,1].unique(), range(len(summary.iloc[:,1].unique())))) | \
            dict(zip(summary.iloc[:,3].unique(), range(len(summary.iloc[:,3].unique()))))
        mapping = lambda x: x.map(map_order)
    else: mapping = None
    
    for i in range(len(group1)):
        for j in range(len(group2)):
            g1 = group1[i]; g2 = group2[j]
            tmp = summary.loc[(summary.iloc[:,0] == g1) & (summary.iloc[:,2] == g2),:]
            tmp = tmp.sort_values(by = [tmp.columns[1], tmp.columns[3]], key = mapping) # important for alignment
            sns.scatterplot(
                tmp,
                x = summary.columns[3], y = summary.columns[1],
                hue = summary.columns[4], hue_norm = (-absmax, absmax),
                palette = 'redblue',
                size = 'Significance', sizes = sizes,
                edgecolor = '.7',
                legend = False,
                ax = ax[i,j]
                )
            
            # dark borders for most significant category
            if tmp.loc[tmp.Significance == sig_label,:].size > 0:
                sns.scatterplot(
                    tmp.loc[tmp.Significance == sig_label,:],
                    x = summary.columns[3], y = summary.columns[1],
                    hue = summary.columns[4], hue_norm = (-absmax, absmax),
                    palette = 'redblue',
                    size = 'Significance', sizes = sizes,
                    edgecolor = 'k', linewidths = 1.5,
                    legend = False,
                    ax = ax[i,j]
                    )
            
            # remove x labels except for last row
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                if i > 0:
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
            
            # diagonal line
            if i == j and autocor: ax[i,j].axline((0,0), slope = 1, c = 'k', zorder = 0)
    
    if autocor: ax[0,0].annotate(annot,(0,0),xytext=(-4,-3.5), rotation = 45)
    
    # colour bar
    norm = mpl.colors.Normalize(vmin = -absmax, vmax = absmax)
    figsize = fig.get_size_inches()
    right_pos = ax[-1,-1].get_position().x1
    cax = fig.add_axes((right_pos+1/figsize[0], 0.4, 0.3/figsize[0], 0.35))
    cax.set_title(summary.columns[4])
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap='redblue'), cax = cax)
    
    # move legend location
    if legend_param:
        handles = [mpl.lines.Line2D([],[], marker = 'o', linewidth = 0, markersize = y**0.5, # Seaborn point size correspond to the square of diameter 
            markeredgecolor = '.7', markerfacecolor = '.7', label = x) for x, y in list(sizes.items())[:-1]] + [
            mpl.lines.Line2D([],[], marker = 'o', linewidth = 0, markersize = list(sizes.values())[-1]**0.5, 
            markeredgecolor = 'k', markeredgewidth = 1.5, markerfacecolor = '.7', label = list(sizes.keys())[-1])]
        fig.legend(handles=handles, title = 'Significance',
                   frameon = False, loc = 'upper left', bbox_to_anchor=(right_pos+0.3/figsize[0], 0.4))
        
    return fig