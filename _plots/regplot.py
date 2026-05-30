#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-11-13

A plotting tool to plot regression plots
'''

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import median_abs_deviation as mad
from .aes import discrete_palette_n, redblue

def _disjoint_groups(df):
    '''
    First column = categorical variable to be regressed separately
    Second column = x axis values to identify discontinuous segments
    greedily group categories into discontinuous groups based on the second column
    '''
    if df.iloc[:,0].unique().size <= 10:
        return [[x] for x in df.iloc[:,0].unique()] # too few categories, can use distinct colours
    else:
        df = df.sort_values(by = df.columns[1]).reset_index(drop = True)
        xlim = df.iloc[:,1].max() - df.iloc[:,1].min()
        out = []
        for cat in df.iloc[:,0].unique(): # groups will be ordered by the second column
            sub_df = df.loc[df.iloc[:,0] == cat, :].reset_index(drop = True)
            cat_time_range = (sub_df.iloc[:,1].min(), sub_df.iloc[:,1].max())
            for i in range(len(out)):
                _, group_ranges = out[i]
                if cat_time_range[0] <= group_ranges[-1][1] + xlim * 0.05: continue
                out[i][0].append(cat)
                out[i][1].append(cat_time_range)
                break
            else:
                out.append([[cat], [cat_time_range]])
        out = [x[0] for x in out]
        return out
    
def _generate_colour_palette(groups):
    '''groups should be a list of lists, generated from _disjoint_groups'''
    n_groups = len(groups)
    palette = discrete_palette_n(n_groups)
    return {cat: palette[i] for i in range(n_groups) for cat in groups[i]}

def _add_regression_axis(
    fig, position, 
    df, x, y, hue, palette, *, 
    df_clipped = None, xlabel = '', ylabel = '', xlabel_groups = True,
    s = 0.5, alpha = 0.2, order = 1, clip_ylim = 0):

    if df_clipped is None: df_clipped = df

    ax = fig.add_axes(position)
    xrange = df[x].max() - df[x].min()
    xlim = (df[x].min() - xrange * 0.05, df[x].max() + xrange * 0.05)
    ax.set_xlim(xlim)
    ax.spines['right'].set_visible(False)
    if xlabel_groups:
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis = 'x', which = 'both', bottom = False, top = True, labelbottom = False, labeltop = True)
        ax.xaxis.set_label_position('top')
    else: ax.spines['top'].set_visible(False)
    if xlabel == '': ax.spines[['top','bottom']].set_visible(False); ax.set_xticks([])
    ax.set_ylim(df[y].quantile(clip_ylim), df[y].quantile(1 - clip_ylim))
    for cat, cat_df in df_clipped.groupby(hue, observed = True):
        sns.regplot(
            data = cat_df,
            x = x, y = y, order = order,
            scatter_kws = dict(s = s, alpha = alpha, color = palette[cat], edgecolor = 'none', rasterized = True),
            line_kws = dict(color = palette[cat], linewidth = 2),
            ax = ax
        )
    for cat, cat_df in df.loc[df.index.difference(df_clipped.index)].groupby(hue, observed = True):
        sns.scatterplot(cat_df, x = x, y = y, color = palette[cat], s = s, alpha = alpha / 2, edgecolor = 'none', ax = ax, rasterized = True)
    ax.set_xlabel(xlabel, fontsize = 12)
    ax.set_ylabel(ylabel, fontsize = 12)
    return ax

def temporal_regplot(
    df,
    x, *y, hue = None, clip_tail = 0.025, height = 3, width = 7,
    xlabel = False, ylabel = None, xlabel_groups = True,
    s = 0.5, alpha = 0.2, clip_ylim = 0,
    order = 1, annotate_corr = True
):
    '''
    df: pd.DataFrame
    x: usually represents time or pseudotime
    y: usually a temporally changing variable
    hue: categorical variable, different categories will be regressed separately
    clip_tail: for each category in hue, clip the top and bottom tail of x axis values by this proportion to avoid outlier effects
        values > 1 are considered MADs and values < 1 are considered quantiles
    xlabel, ylabel: axis labels, specify False to hide, None to use the column name
    if xlabel is False, the horizontal spine will be hidden
    xlabel_groups: whether to show x axis labels for each group separately, only applicable when hue is not None
    s: point size
    alpha: point transparency
    clip_ylim: proportion of y-axis limits to clip at the top and bottom
    '''

    # region: argument sense checks and basic config
    df = df.dropna()
    df[x] = df[x].astype(float) # ensure x is numeric
    if isinstance(ylabel, str): ylabel = [ylabel] # allow y to be a single column name
    y = list(y)
    if xlabel is None: xlabel = x.replace('_',' ')
    if ylabel is None: ylabel = [col.replace('_',' ') for col in y]
    if xlabel is False: xlabel = ''
    if ylabel is False: ylabel = ['' for _ in y]
    if not isinstance(xlabel, str): raise ValueError('xlabel should be a string, None, or False')
    if not isinstance(ylabel, list): raise ValueError('ylabel should be a string, a list of strings, None, or False')
    assert len(y) == len(ylabel), 'length of y and ylabel should match'
    sns.set_style('ticks')
    if order == 2: annotate_corr = False # disable correlation annotation for quadratic regression
    if hue is None: xlabel_groups = False # no groups if hue is None

    total_height = height * len(y)
    xrange = df[x].max() - df[x].min()
    xlim = (df[x].min() - xrange * 0.05, df[x].max() + xrange * 0.05)
    y = y[::-1]; ylabel = ylabel[::-1] # reverse the order of y and ylabel so that the first one is plotted at the top
    # endregion

    # region: simple regression plots
    if hue is None:
        fig = plt.Figure(figsize = (width, total_height+0.3))
        current_height = 0.3
        for col, col_label in zip(y, ylabel):
            ax_position = (0.3/width, current_height/total_height, (width-0.6)/width, (height-0.3)/total_height)
            ax = fig.add_axes(ax_position)
            ax.set_xlim(xlim)
            ax.spines[['top','right']].set_visible(False)
            if xlabel == '': ax.spines['bottom'].set_visible(False); ax.set_xticks([])
            ax.set_ylim(df[col].quantile(clip_ylim), df[col].quantile(1 - clip_ylim))

            sns.regplot(
                data = df,
                x = x, y = col, order = order, 
                scatter_kws = dict(s = s, alpha = alpha, color = '0.8', edgecolor = 'none', rasterized = True),
                line_kws = dict(color = '0.8', linewidth = 2),
                ax = ax
            )
            ax.set_xlabel(xlabel, fontsize = 12)
            ax.set_ylabel(ylabel, fontsize = 12)

            # annotate correlation value
            if annotate_corr:
                r = df[x].corr(df[col])
                ax.text(ax_position.x0+0.5/width, ax_position.y1-0.5/height,
                        f'r = {r:.3f}', fontsize = 12, ha = 'left', va = 'top', color = 'k')
        return fig
    # endregion

    # clip tails 
    if clip_tail > 0:
        df_list = []
        for cat, sub_df in df.groupby(hue, observed = True):
            if clip_tail > 1:
                lower = sub_df[x].median() - clip_tail * mad(sub_df[x])
                upper = sub_df[x].median() + clip_tail * mad(sub_df[x])
            else:
                lower = sub_df[x].quantile(clip_tail)
                upper = sub_df[x].quantile(1 - clip_tail)
            df_list.append(sub_df.loc[(sub_df[x] >= lower) & (sub_df[x] <= upper), :])
        df_clipped = pd.concat(df_list, axis = 0)
    else: df_clipped = df.copy()

    # separate the categories into disjoint groups
    groups = _disjoint_groups(df_clipped[[hue, x]])
    palette = _generate_colour_palette(groups)

    # if xlabels are needed for each group, then we need another axis for annotations
    if xlabel_groups:
        annot_height = len(groups)/4 # 0.25 inches per group
        total_height += annot_height
        fig = plt.Figure(figsize = (width, total_height)) # 0.25 inches per group
        # specify details for the by-group annotation
        annot_axis = fig.add_axes((
            0.3/width, 0.05/total_height, 
            (width-0.6)/width, (annot_height-0.05)/total_height
        ))
        annot_axis.spines[['top','right','bottom','left']].set_visible(False)
        annot_axis.set_xticks([]); annot_axis.set_yticks([])
        annot_height = annot_height # 0.33 inches per group
        annot_axis.set_ylim(-annot_height-0.1, 0); annot_axis.set_xlim(xlim)
        current_height = -0.05
        for group in groups:
            for cat in group:
                cat_clipped = df_clipped.loc[df_clipped[hue] == cat, :]
                annot_axis.hlines(current_height, color = palette[cat], linewidth = 2,
                    xmax = cat_clipped[x].max(), xmin = cat_clipped[x].min())
                if annotate_corr and len(y) == 1: # only one correlation value can be annotated
                    cat_full = df.loc[df[hue] == cat, :]
                    r = cat_full[x].corr(cat_full[y[0]])
                    annot_axis.text((cat_clipped[x].min() + cat_clipped[x].max()) / 2, current_height - 0.05,
                        f'{cat}: r={r:.3f}', fontsize = 12, ha = 'center', va = 'top', color = palette[cat])
                else: annot_axis.text((cat_clipped[x].min() + cat_clipped[x].max()) / 2, current_height - 0.05,
                    cat, fontsize = 12, ha = 'center', va = 'top', color = palette[cat])
            current_height -= 1/4
    else: fig = plt.Figure(figsize = (width, total_height))

    # add regression axes
    current_height = annot_height if xlabel_groups else 0
    for col, col_label in zip(y, ylabel):
        ax_position = (0.3/width, current_height/total_height, (width-0.6)/width, (height-0.3)/total_height)
        ax_xlabel = xlabel if col == y[-1] else '' # only show xlabel for the top plot
        _add_regression_axis(
            fig, ax_position, df, x, col, hue, palette,
            df_clipped = df_clipped, xlabel = ax_xlabel, ylabel = col_label,
            xlabel_groups = xlabel_groups, s = s, alpha = alpha, order = order, clip_ylim = clip_ylim
        )
        current_height += height

    return fig
    