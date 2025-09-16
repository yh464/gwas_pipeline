#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-09-16

Plotting aesthetics
'''

import matplotlib as mpl

# colour bars
redblue = mpl.colors.LinearSegmentedColormap(
  'redblue',
  dict(red = ((0,0,0),(1/2,1,1),(1,.8,.8)),
      green = ((0,0,0),(1/2,1,1),(1,0,0)),
      blue = ((0,.8,.8),(1/2,1,1),(1,0,0))),
  1024
  )