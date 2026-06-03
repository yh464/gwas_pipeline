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
from .aes import register_palettes, discrete_palette

def _add_rep_axis(fig, rep = 'UMAP'):
  # UMAP axes
  figsize = fig.get_size_inches()
  repax = fig.add_axes((0.25/figsize[0], 0.25/figsize[1], 0.5/figsize[0], 0.5/figsize[1]))
  repax.spines[["top", "right"]].set_visible(False)
  repax.set_xticks([]); repax.set_yticks([])
  repax.set_xlabel(f'{rep}1', fontsize = 8); repax.set_ylabel(f'{rep}2', fontsize = 8)
  return fig

def scatterplot_noaxis(x, y, v, *, full = True, palette = None, s = 'auto', rep = 'UMAP', vname = '', vmin = None, vmax = None, **kwargs):
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
  if not full and df.shape[0] > 500000: df = df.sample(500000, random_state = 0) # for large datasets, sample points to speed up plotting

  # colour palette
  if palette != None:
    register_palettes(palette)
    use_palette = palette
  else:
    if v.dtype.name == 'category' or v.dtype == object:
      use_palette = discrete_palette(v.unique())
    else:
      if np.nanmax(v) <= 0: from .aes import greyblue_alpha0; use_palette = greyblue_alpha0.name; register_palettes(greyblue_alpha0)
      elif np.nanmin(v) >= 0: from .aes import redgrey_alpha0; use_palette = redgrey_alpha0.name; register_palettes(redgrey_alpha0)
      else: from .aes import redblue_alpha; register_palettes(redblue_alpha); use_palette = redblue_alpha.name
  legend = 'auto' if (v.dtype.name == 'category' or v.dtype == object) else False
  if not (v.dtype.name == 'category' or v.dtype == object):
    if np.nanmax(v) <= 0: vmin = max(np.nanquantile(v, 0.05)*1.5,np.nanmin(v)) if vmin == None else vmin; vmax = 0 if vmax == None else vmax
    elif np.nanmin(v) >= 0: vmin = 0 if vmin == None else vmin; vmax = min(np.nanquantile(v,0.95),np.nanmax(v)) if vmax == None else vmax
    else:
      absv = np.abs(v)
      vmin = -max(np.nanquantile(absv, 0.95)*1.5,np.nanmax(absv)) if vmin == None else vmin
      vmax = max(np.nanquantile(absv, 0.95)*1.5,np.nanmax(absv)) if vmax == None else vmax
      vmax, vmin = max(abs(vmin), abs(vmax)), -max(abs(vmin), abs(vmax)) # ensure vmin and vmax are symmetric around 0
    kwargs['hue_norm'] = mpl.colors.Normalize(vmin = vmin, vmax = vmax)

  # main plot
  fig = plt.figure(figsize = (5.3,5))
  ax = fig.add_axes((0.3/5.3, 0.3/5, 4.4/5.3, 4.4/5))
  if s == 'auto':
    # estimate distance between points in the x-y plane
    from scipy.spatial.distance import pdist
    if df.shape[0] > 1 and df.shape[0] < 50000:
      dists = pdist(df[['x','y']].values)
      min_dist = np.nanquantile(dists, 0.05)
      axis_range = max(df['x'].max() - df['x'].min(), df['y'].max() - df['y'].min())
      s = (min_dist/axis_range)**2 * 144
    elif df.shape[0] >= 50000: s = 0.1
    else: s = 10
    s = min(max(s, 0.1), 36) # set a reasonable range for point size
  sns.scatterplot(data = df, x = 'x', y = 'y', hue = 'v', palette = use_palette, s = s, ax = ax, edgecolor = None, linewidth = 0, 
      legend = legend, rasterized = True, **kwargs)
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
  
  if isinstance(rep, str): fig = _add_rep_axis(fig, rep.upper())
  return fig

def scatterplot_adata(adata, v, rep = 'umap', **kwargs):
  '''
  Scatterplot from anndata object
  Input:
    adata: anndata object with obsm[rep] containing coordinates
    v: variable name in adata.obs or adata.var to colour by
    rep: representation in adata.obsm to use for coordinates
    all **kwargs are passed to scatterplot_noaxis()
  '''
  rep = rep[2:] if rep.startswith('X_') else rep
  if rep.lower() == 'umap':
    coord_key = 'X_umap'; rep_axis = 'UMAP'
  elif rep.lower() == 'pca':
    coord_key = 'X_pca'; rep_axis = 'PC'
  elif rep.lower() == 'tsne':
    coord_key = 'X_tsne'; rep_axis = 'tSNE'
  elif rep.lower() == 'spatial':
    coord_key = 'spatial'; rep_axis = False
  elif rep in adata.obsm.keys():
    coord_key = rep; rep_axis = rep.upper().replace('X_','')
  else:
    raise ValueError(f'Reduced-dimension representation {rep} not in the anndata object')
  x = adata.obsm[coord_key][:,0]
  y = adata.obsm[coord_key][:,1]
  
  if isinstance(v, str) and v in adata.obs.columns.tolist():
    val = adata.obs[v]; vname = v
  elif isinstance(v, str) and v in adata.var.index.tolist():
    gene_idx = adata.var.index.get_loc(v)
    val = adata.X[:,gene_idx].toarray().flatten() if hasattr(adata.X, 'toarray') else adata.X[:,gene_idx].flatten()
    vname = v
  else:
    val = v; vname = ''
  
  fig = scatterplot_noaxis(x, y, val, rep = rep_axis, vname = vname, **kwargs)
  return fig