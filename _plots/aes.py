#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-09-16

Plotting aesthetics
'''

import matplotlib as mpl
import seaborn as sns
import math

# colour palettes
def register_palettes(*palettes):
  for p in palettes:
    try: mpl.colormaps.register(p)
    except: pass

redblue = mpl.colors.LinearSegmentedColormap(
  'redblue',
  dict(red = ((0,0,0),(1/2,1,1),(1,.8,.8)),
      green = ((0,0,0),(1/2,1,1),(1,0,0)),
      blue = ((0,.8,.8),(1/2,1,1),(1,0,0))),
  1024
)

redblue_alpha = mpl.colors.LinearSegmentedColormap(
  'redblue_alpha',
  dict(red = ((0,0,0),(1/2,1,1),(1,.8,.8)),
      green = ((0,0,0),(1/2,1,1),(1,0,0)),
      blue = ((0,.8,.8),(1/2,1,1),(1,0,0)),
      alpha = ((0,1,1),(1/2,0.2,0.2),(1,1,1))),
  1024
)

redwhite = mpl.colors.LinearSegmentedColormap(
  'redwhite',
  dict(red = ((0,1,1),(1,0.8,0.8)),
       green = ((0,1,1),(1,0,0)),
       blue = ((0,1,1),(1,0,0))),
  1024
)

redgrey = mpl.colors.LinearSegmentedColormap(
  'redgrey',
  dict(red = ((0,0.8,0.8),(1,0.8,0.8)),
       green = ((0,0.8,0.8),(1,0,0)),
       blue = ((0,0.8,0.8),(1,0,0))),
  1024
  )

redgrey_alpha = mpl.colors.LinearSegmentedColormap(
  'redgrey_alpha',
  dict(red = ((0,0.8,0.8),(1,0.8,0.8)),
       green = ((0,0.8,0.8),(1,0,0)),
       blue = ((0,0.8,0.8),(1,0,0)),
       alpha = ((0,0.2,0.2),(1,1,1))),
  1024
  )

whiteblue = mpl.colors.LinearSegmentedColormap(
  'whiteblue',
  dict(red = ((0,0,0),(1,1,1)),
       green = ((0,0,0),(1,1,1)),
       blue = ((0,0.8,0.8),(1,1,1))),
  1024,
)

greyblue = mpl.colors.LinearSegmentedColormap(
  'greyblue',
  dict(red = ((0,0,0),(1,.8,.8)),
       green = ((0,0,0),(1,.8,.8)),
       blue = ((0,0.8,0.8),(1,.8,.8))),
  1024,
)

greyblue_alpha = mpl.colors.LinearSegmentedColormap(
  'greyblue_alpha',
  dict(red = ((0,0,0),(1,.8,.8)),
       green = ((0,0,0),(1,.8,.8)),
       blue = ((0,0.8,0.8),(1,.8,.8)),
       alpha = ((0,1,1),(1,.2,.2))),
  1024,
)

def discrete_palette(n):
  '''Default colour palette to use for discrete mapping'''
  if n <= 10: return sns.color_palette('muted', n)
  else:
    palette = sns.color_palette('husl', n) # reorder to maximise distance between adjacent colours
    split_half = math.ceil(n / 2)
    out = []
    for i in range(split_half):
      out.append(palette[i])
      if i + split_half < n:
        out.append(palette[i + split_half])
    return out