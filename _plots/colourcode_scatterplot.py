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
from .aes import register_palettes

def _add_rep_axis(fig, rep = 'UMAP'):
  # UMAP axes
  figsize = fig.get_size_inches()
  repax = fig.add_axes((0.25/figsize[0], 0.25/figsize[1], 0.5/figsize[0], 0.5/figsize[1]))
  repax.spines[["top", "right"]].set_visible(False)
  repax.set_xticks([]); repax.set_yticks([])
  repax.set_xlabel(f'{rep}1', fontsize = 8); repax.set_ylabel(f'{rep}2', fontsize = 8)
  return fig

def scatterplot_noaxis(x, y, v, palette = None, s = 0.1, rep = 'UMAP', vname = '', vmin = None, vmax = None, **kwargs):
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
      if np.nanmax(v) <= 0: from .aes import whiteblue; use_palette = whiteblue.name; register_palettes(whiteblue)
      elif np.nanmin(v) >= 0: from .aes import redgrey_alpha; use_palette = redgrey_alpha.name; register_palettes(redgrey_alpha)
      else: from .aes import redblue_alpha; register_palettes(redblue_alpha); use_palette = redblue_alpha.name
  legend = 'auto' if (v.dtype.name == 'category' or v.dtype == object) else False
  if not (v.dtype.name == 'category' or v.dtype == object):
    if np.nanmax(v) <= 0: vmin = max(np.nanquantile(v, 0.05)*1.5,np.nanmin(v)) if vmin == None else vmin; vmax = 0 if vmax == None else vmax
    elif np.nanmin(v) >= 0: vmin = 0 if vmin == None else vmin; vmax = min(np.nanquantile(v,0.95),np.nanmax(v)) if vmax == None else vmax
    else:
      absv = np.abs(v) 
      vmin = -max(np.nanquantile(absv, 0.95)*1.5,np.nanmax(absv)) if vmin == None else vmin
      vmax = max(np.nanquantile(absv, 0.95)*1.5,np.nanmax(absv)) if vmax == None else vmax
    kwargs['hue_norm'] = mpl.colors.Normalize(vmin = vmin, vmax = vmax)

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
    sns.move_legend(ax, "center left", bbox_to_anchor=(1, 0.5), title = 'Group', frameon = False, ncols = np.ceil(len(v.unique())/20))
    # increase point size in legend
    handles = ax.get_legend().legend_handles
    for handle in handles: handle.set_markersize(5)
  
  if rep != False: fig = _add_rep_axis(fig, rep)
  return fig
