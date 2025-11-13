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
from .aes import discrete_palette, redblue

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
    palette = discrete_palette(n_groups)
    return {cat: palette[i] for i in range(n_groups) for cat in groups[i]}

def temporal_regplot(
    df,
    x, y, hue = None, clip_tail = 0.025,
    xlabel = None, ylabel = None, xlabel_groups = True,
    s = 0.5, alpha = 0.2, 
    order = 1, annotate_corr = True
):
    '''
    df: pd.DataFrame
    x: usually represents time or pseudotime
    y: usually a temporally changing variable
    hue: categorical variable, different categories will be regressed separately
    clip_tail: for each category in hue, clip the top and bottom tail of x axis values by this proportion to avoid outlier effects
        values > 1 are considered MADs and values < 1 are considered quantiles
    xlabel, ylabel: axis labels
    xlabel_groups: whether to show x axis labels for each group separately
    s: point size
    alpha: point transparency
    '''

    df = df.dropna()
    if xlabel is None: xlabel = x.replace('_',' ')
    if ylabel is None: ylabel = y.replace('_',' ')
    sns.set_style('ticks')
    if order == 2: annotate_corr = False # disable correlation annotation for quadratic regression

    # simple regression plots
    if hue is None:
        fig, ax = plt.subplots(figsize = (7, 4))
        ax_position = ax.get_position()
        ax.spines[['top','right','bottom']].set_visible(False)
        ax.set_xticks([])
        ax.set_ylim(df[y].min(), df[y].max())
        sns.regplot(
            data = df,
            x = x, y = y, order = order, 
            scatter_kws = dict(s = s, alpha = alpha, color = '0.8', edgecolor = 'none'),
            line_kws = dict(color = '0.8', linewidth = 2),
            ax = ax
        )
        ax.set_xlabel(xlabel, fontsize = 12)
        ax.set_ylabel(ylabel, fontsize = 12)

        # annotate correlation value
        if annotate_corr:
            r = df[x].corr(df[y])
            ax.text(ax_position.x0+0.5/figsize[0], ax_position.y1-0.5/figsize[1],
                    f'r = {r:.3f}', fontsize = 12, ha = 'left', va = 'top', color = 'k')
        return fig

    # clip tails 
    if clip_tail > 0:
        df_list = []
        excluded_list = []
        for cat, sub_df in df.groupby(hue, observed = True):
            if clip_tail > 1:
                lower = sub_df[x].median() - clip_tail * mad(sub_df[x])
                upper = sub_df[x].median() + clip_tail * mad(sub_df[x])
            else:
                lower = sub_df[x].quantile(clip_tail)
                upper = sub_df[x].quantile(1 - clip_tail)
            df_list.append(sub_df.loc[(sub_df[x] >= lower) & (sub_df[x] <= upper), :])
            excluded_list.append(sub_df.loc[(sub_df[x] < lower) | (sub_df[x] > upper), :])
        df_clipped = pd.concat(df_list, axis = 0)
        df_excluded = pd.concat(excluded_list, axis = 0)
    else:
        df_clipped = df.copy()
        df_excluded = pd.DataFrame(columns = df.columns, index = [])

    # separate the categories into disjoint groups
    groups = _disjoint_groups(df_clipped[[hue, x]])
    palette = _generate_colour_palette(groups)


    # if xlabels are needed for each group, then we need another axis for annotations
    if xlabel_groups:
        annot_height = len(groups)/4 # 0.33 inches per group
        fig = plt.Figure(figsize = (7, 4 + annot_height)) # 0.33 inches per group
        ax = fig.add_axes((0.3/7, annot_height/ (4 + annot_height), 6.4/7, 3.7/(4 + annot_height)))
        annot_axis = fig.add_axes((0.3/7, 0.05/ (4 + annot_height), 6.4/7, (annot_height-0.05)/(4 + annot_height)))
        annot_axis.spines[['top','right','bottom','left']].set_visible(False)
        annot_axis.set_xticks([]); annot_axis.set_yticks([])
    else:
        fig = plt.Figure(figsize = (7, 4))
        ax = fig.add_axes((0.3/7, 0.3/4, 6.4/7, 3.4/4))
    ax_position = ax.get_position()
    figsize = fig.get_size_inches()
    ax.spines[['top','right','bottom']].set_visible(False)
    ax.set_xticks([])
    ax.set_ylim(df[y].min(), df[y].max())

    # plot regression lines for each category with tails clipped
    for cat, cat_df in df_clipped.groupby(hue, observed = True):
        sns.regplot(
            data = cat_df,
            x = x, y = y, order = order,
            scatter_kws = dict(s = s, alpha = alpha, color = palette[cat], edgecolor = 'none'),
            line_kws = dict(color = palette[cat], linewidth = 2),
            ax = ax
        )
    # plot clipped tails with the same colour but more transparent
    if df_excluded.shape[0] > 0:
        for cat, cat_df in df_excluded.groupby(hue, observed = True):
            sns.scatterplot(cat_df, x = x, y = y, color = palette[cat], s = s, alpha = alpha / 2, edgecolor = 'none', ax = ax)

    # annotate correlation values for each category
    if xlabel_groups:
        annot_height = annot_height # 0.33 inches per group
        annot_axis.set_ylim(-annot_height-0.1, 0); annot_axis.set_xlim(ax.get_xlim())
        current_height = -0.05
        for group in groups:
            for cat in group:
                cat_clipped = df_clipped.loc[df_clipped[hue] == cat, :]
                annot_axis.hlines(current_height, color = palette[cat], linewidth = 2,
                    xmax = cat_clipped[x].max(), xmin = cat_clipped[x].min())
                if annotate_corr:
                    cat_full = df.loc[df[hue] == cat, :]
                    r = cat_full[x].corr(cat_full[y])
                    annot_axis.text((cat_clipped[x].min() + cat_clipped[x].max()) / 2, current_height - 0.05,
                        f'{cat}: r={r:.3f}', fontsize = 12, ha = 'center', va = 'top', color = palette[cat])
                else: annot_axis.text((cat_clipped[x].min() + cat_clipped[x].max()) / 2, current_height - 0.05,
                    cat, fontsize = 12, ha = 'center', va = 'top', color = palette[cat])
            current_height -= 1/4
        ax.set_xlabel('')
    else: ax.set_xlabel(xlabel, fontsize = 12)
    ax.set_ylabel(ylabel, fontsize = 12)
    return fig
    