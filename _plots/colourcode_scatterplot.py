#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2025-09-16

A flexible framework to plot scatterplots where coordinates are (spatial, UMAP, etc) representations
and colour codes represent the property of interest
'''

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import scanpy as sc
from .aes import register_palettes

def _add_rep_axis(fig, rep = 'UMAP'):
  # UMAP axes
  figsize = fig.get_size_inches()
  repax = fig.add_axes((0.25/figsize[0], 0.25/figsize[1], 0.5/figsize[0], 0.5/figsize[1]))
  repax.spines[["top", "right"]].set_visible(False)
  repax.set_xticks([]); repax.set_yticks([])
  repax.set_xlabel(f'{rep}1', fontsize = 8); repax.set_ylabel(f'{rep}2', fontsize = 8)
  return fig

def scatterplot_noaxis(x, y, v, palette = None, s = 0.1, rep = 'UMAP', vname = '', **kwargs):
  '''
  Scatterplot without axes
  Input:
    x, y: coordinates
    v: values for colouring (continuous or categorical)
    palette: a matplotlib.colors.Colormap
    s: point size
    all **kwargs are passed to sns.scatterplot()
  '''
  df = pd.DataFrame(dict(x = x, y = y, v = v)).dropna()

  # colour palette
  if palette != None:
    register_palettes(palette)
    use_palette = palette.name
  else:
    if v.dtype.name == 'category' or v.dtype == object:
      use_palette = sns.color_palette('husl', n_colors = len(v.unique()))
    else:
      from .aes import redblue
      register_palettes(redblue)
      use_palette = redblue.name
  legend = 'auto' if (v.dtype.name == 'category' or v.dtype == object) else False
  if not (v.dtype.name == 'category' or v.dtype == object):
    if (v < 0).all(): kwargs['hue_norm'] = mpl.colors.Normalize(vmin = np.min(v), vmax = 0)
    elif (v > 0).all(): kwargs['hue_norm'] = mpl.colors.Normalize(vmin = 0, vmax = np.max(v))
    else: kwargs['hue_norm'] = mpl.colors.Normalize(vmin = -np.max(np.abs(v)), vmax = np.max(np.abs(v)))

  # main plot
  fig = plt.figure(figsize = (5.3,5))
  ax = fig.add_axes((0.3/5.3, 0.3/5, 4.4/5.3, 4.4/5))
  sns.scatterplot(data = df, x = 'x', y = 'y', hue = 'v', palette = use_palette, s = s, ax = ax, edgecolor = None, linewidth = 0, legend = legend, **kwargs)
  ax.axis('off')
  ax.set_aspect('equal')
  
  # colour bar
  if not (v.dtype.name == 'category' or v.dtype == object):
    norm = kwargs['hue_norm']
    cax = fig.add_axes((4.8/5.3, 0.3/5, 0.2/5.3, 4.4/5))
    cax.set_title(vname)
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=use_palette), cax = cax)
  # legend
  else:
    sns.move_legend(ax, "center left", bbox_to_anchor=(1, 0.5), title = 'Group', frameon = False)
  
  if rep != False: fig = _add_rep_axis(fig, rep)
  return fig
