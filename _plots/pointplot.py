#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2026-05-23

A generalised plotting tool to plot pointplots with error bars
'''

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from .tidy import capitalise, get_fdr
from .aes import discrete_palette

def summary_pointplot(summary, xgroup = None, x = None, hue = None, y = None, se = None, 
    sort = False, p_threshold: list[float] = [], sig_col = None, xlabel = True, ylabel = None):

    '''
    Default input format: long format pd.DataFrame, compatible with corr_heatmap
        1st column: ignored
        2nd column: individual label (hue), usually a phenotype
        3rd column: group label (x axis, as xlabel)
        4th column: individual label (x axis, as xticklabels), usually a cell type
        5th column: statistic value (y axis)
        6th column: standard error (for error bars) (prioritise 'se' column if exists)
        and may contain a 'p' and 'q' column somewhere in the data frame
    '''

    # input check
    if xgroup is None: xgroup = summary.columns[2]
    if x is None: x = summary.columns[3]
    if hue is None: hue = summary.columns[1]
    if y is None: y = summary.columns[4]
    if se is None: se = 'se' if 'se' in summary.columns else summary.columns[5]
    assert not any([tmp is None for tmp in [xgroup, x, hue, y, se]]), 'xgroup, x, hue, y and se must be specified or automatically inferred from the first 6 columns of the input data frame'
    if sort: summary = summary.sort_values([xgroup, x, hue])
    summary, sig_label = get_fdr(summary, group_by = [], sig_col = sig_col, p_threshold = p_threshold)
    summary[xgroup] = capitalise(summary[xgroup])
    summary[x] = capitalise(summary[x])
    summary[hue] = capitalise(summary[hue])

    # style sheet
    sns.set_theme(style = 'ticks', rc = {
        'axes.spines.right': False, 'axes.spines.top': False, 'axes.spines.bottom': False,
        'xtick.bottom': False
        })
    plt.tick_params(axis = 'x', rotation = 90)

    # determine figure size
    groups = summary[xgroup].unique()
    counts = [summary.loc[summary[xgroup] == xgroup_i, x].unique().size for xgroup_i in groups]
    hues = summary[hue].unique()
    fig, ax = plt.subplots(1, len(groups), width_ratios = counts, sharey = False, squeeze = True,
        figsize = (sum(counts) * 0.3, 3), gridspec_kw = {'wspace': 0.1})
    if len(groups) == 1: ax = np.array([ax]) # ensure ax is always an array for consistent indexing
    ylim = (min(summary[y] - summary[se] * 1.96), max(summary[y] + summary[se] * 1.96))
    
    # aesthetics
    palette = discrete_palette(hues) # need to use dict mapping as sig and non-sig are plotted separately
    if not sort:
        map_order = dict(zip(summary[hue].unique(), range(len(summary[hue].unique())))) | \
            dict(zip(summary[x].unique(), range(len(summary[x].unique()))))
        mapping = lambda z: z.map(map_order)
    else: mapping = None

    for i, group in enumerate(groups):
        tmp = summary.loc[summary[xgroup] == group,:]
        tmp = tmp.sort_values(by = [x, hue], key = mapping) # important for alignment

        # manually specify x coordinates for each point
        dodge_width = 0.05 * tmp[hue].unique().size
        col_width = max(1, dodge_width * 2 + 0.2)
        tmp['x_pos'] = tmp[x].map(dict(zip(tmp[x].unique(), range(tmp[x].unique().size)))) * col_width
        tmp['x_pos'] += tmp[hue].map(dict(zip(tmp[hue].unique(), np.linspace(-dodge_width, dodge_width, tmp[hue].unique().size)))) # dodge

        # use scatterplot for marker size control
        tmp['Sig_bin'] = tmp['Significance'] == sig_label
        sns.scatterplot(tmp, x = 'x_pos', y = y, hue = hue, palette = palette, ax = ax[i],
            size = tmp['Sig_bin'].map({True: 250, False: 25}), marker = 'o', edgecolors = 'none', legend = False)

        for _, row in tmp.iterrows():
            ax[i].errorbar(
                row['x_pos'], row[y], yerr = row[se] * 1.96, fmt = 'none', 
                ecolor = palette[row[hue]], 
                elinewidth = 3 if row['Significance'] == sig_label else 1,
                capsize = 0
            )
        
        ax[i].set_xlim(tmp['x_pos'].min() - 0.5*col_width, tmp['x_pos'].max() + 0.5*col_width)
        ax[i].set_xlabel(group if xlabel else '', fontsize = 12)
        ax[i].set_ylim(ylim)
        if i > 0: 
            ax[i].spines['left'].set_visible(False)
            ax[i].set_ylabel('')
            ax[i].set_yticks([])
        else:
            ax[i].set_ylabel(y if ylabel is None else ylabel, fontsize = 12)
            ax[i].spines['left'].set(position = ('outward', 10))
        ax[i].set_xticks(np.array(range(tmp[x].unique().size)) * col_width, tmp[x].unique(), rotation = 90)
        ax[i].axhline(0, color = '0.7', linestyle = '--', linewidth = 1.5, zorder = 0)
    
    # legend
    # colour code for all groups specified in the hue variable
    handles = [mpl.lines.Line2D([],[], marker = 'o', linestyle = 'none', markersize = 5, 
        markerfacecolor = col, markeredgecolor = col, label = cat) for cat, col in palette.items()]
    leg1 = ax[-1].legend(handles = handles, title = hue, loc = 'lower left', frameon = False, bbox_to_anchor = (1, 0.4))

    # size, linestyle and linewidth code for significant groups, only for the most significant category
    # legend should show a marker and a line
    if sig_label != '':
        handles = [
            mpl.lines.Line2D([0, 22], [3.85, 3.85], color = 'k', linewidth = 3, marker = 'o', markersize = 5, 
                label = sig_label, solid_capstyle = 'butt'),
            mpl.lines.Line2D([0, 22], [3.85, 3.85], color = 'k', linewidth = 1, marker = 'o', markersize = np.sqrt(5), 
                label = f'Nominal/NS', solid_capstyle = 'butt')
        ]
        ax[-1].legend(handles = handles, title = 'Significance', loc = 'upper left', frameon = False, bbox_to_anchor = (1, 0.4))
        ax[-1].add_artist(leg1)
    return fig