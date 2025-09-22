#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2025-09-16

Scatterplot of cell-level datasets in UMAP space
'''

from atexit import register
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
  repax = fig.add_axes((0.01, 0.01, 0.5/figsize[0], 0.5/figsize[1]))
  repax.spines[["top", "right"]].set_visible(False)
  repax.set_xticks([]); repax.set_yticks([])
  repax.set_xlabel(f'{rep}1', fontsize = 8); repax.set_ylabel(f'{rep}2', fontsize = 8)
  return fig

def scatterplot_noaxis(x, y, v, palette = None, s = 0.1, rep = 'UMAP', **kwargs):
  '''
  Scatterplot without axes
  Input:
    x, y: coordinates
    v: values for colouring (continuous or categorical)
    palette: a matplotlib.colors.Colormap
    s: point size
    all **kwargs are passed to sns.scatterplot()
  '''
  df = pd.DataFrame(dict(x = x, y = y, v = v))

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
    kwargs['hue_norm'] = mpl.colors.Normalize(vmin = -np.nanmax(np.abs(v)), vmax = np.nanmax(np.abs(v)))

  # main plot
  fig = plt.figure(figsize = (5.3,5))
  ax = fig.add_axes((0.1/5.3, 0.1/5, 4.8/5.3, 4.8/5))
  sns.scatterplot(data = df, x = 'x', y = 'y', hue = 'v', palette = use_palette, s = s, ax = ax, edgecolor = None, linewidth = 0, legend = legend, **kwargs)
  ax.axis('off')
  ax.set_aspect('equal')
  
  # colour bar
  if not (v.dtype.name == 'category' or v.dtype == object):
    norm = mpl.colors.Normalize(vmin = -np.nanmax(np.abs(v)), vmax = np.nanmax(np.abs(v)))
    cax = fig.add_axes((5/5.3, 0.1/5, 0.2/5.3, 4.8/5))
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=use_palette), cax = cax)
  # legend
  else:
    sns.move_legend(ax, "center left", bbox_to_anchor=(1, 0.5), title = 'Group', frameon = False)

  fig = _add_rep_axis(fig, rep)
  return fig
