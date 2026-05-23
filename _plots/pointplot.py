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
    sort = False, p_threshold: list[float] = [], sig_col = None, xlabel = True):

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
    fig, ax = plt.subplots(1, len(groups), width_ratios = counts, sharey = True, squeeze = True,
        figsize = (sum(counts)/3, 3))
    
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

        # initialise x axis positions using a blank plot
        sns.scatterplot(tmp, x = x, y = y, marker = '', ax = ax[i], legend = False)

        # a hack to plot errorbars with just summary data
        tmp_lower = tmp.copy(); tmp_lower[y] = tmp[y] - tmp[se]
        tmp_upper = tmp.copy(); tmp_upper[y] = tmp[y] + tmp[se]
        tmp = pd.concat([tmp, tmp_lower, tmp_upper], axis = 0)

        tmp_sig = tmp.loc[tmp['Significance'] == sig_label, :]
        tmp_ns = tmp.loc[tmp['Significance'] != sig_label, :]

        sns.pointplot(tmp_ns, x = x, y = y, hue = hue, palette = palette, ax = ax[i], dodge = 0.2,
            linestyle = 'none', errorbar = lambda x: (x.min(), x.max()), legend = False,
            markersize = 5, err_kws = {'linewidth': 1, 'capsize': 0, 'linestyle': '--'})
        sns.pointplot(tmp_sig, x = x, y = y, hue = hue, palette = palette, ax = ax[i], dodge = 0.2,
            linestyle = 'none', errorbar = lambda x: (x.min(), x.max()), legend = False,
            markersize = 25, err_kws = {'linewidth': 2, 'capsize': 0})
        
        ax[i].set_xlabel(group if xlabel else '', fontsize = 12)
        if i > 0: 
            ax[i].spines['left'].set_visible(False)
            ax[i].set_ylabel('')
            ax[i].yticks([])
        ax[i].axhline(0, color = '0.7', linestyle = '--', linewidth = 1.5)
    
    # legend
    right_pos = ax[-1].get_position().x1
    figsize = fig.get_size_inches()

    # colour code for all groups specified in the hue variable
    handles = [mpl.lines.Line2D([],[], marker = 'o', linetype = 'none', markersize = 5, 
        markerfacecolor = col, markeredgecolor = col, label = cat) for cat, col in palette.items()]
    fig.legend(handles = handles, title = hue, loc = 'lower left', frameon = False, bbox_to_anchor = (right_pos+0.3/figsize[0], 0.4))

    # size, linestyle and linewidth code for significant groups, only for the most significant category
    # legend should show a marker and a line
    if sig_label != '':
        handles = [
            mpl.lines.Line2D([0, 22], [3.85, 3.85], color = 'k', likestyle = '-', linewidth = 2, label = sig_label),
            mpl.lines.Line2D([0, 22], [3.85, 3.85], color = 'k', likestyle = '--', linewidth = 1, label = f'Nominal/NS')
        ]
        fig.legend(handles = handles, title = 'Significance', loc = 'upper left', frameon = False, bbox_to_anchor = (right_pos+0.3/figsize[0], 0.4))
    