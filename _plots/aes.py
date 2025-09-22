#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-09-16

Plotting aesthetics
'''

import matplotlib as mpl

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

whiteblue = mpl.colors.LinearSegmentedColormap(
  'whiteblue',
  dict(red = ((0,0,0),(1,1,1)),
       green = ((0,0,0),(1,1,1)),
       blue = ((0,0.8,0.8),(1,1,1))),
  1024,
)

