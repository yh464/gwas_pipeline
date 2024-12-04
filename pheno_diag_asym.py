#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2024-11-27

This is a diagnostic script to understand asymmetry phenotypes of the network
Intended to be one-off
'''

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

os.chdir('/rds/user/yh464/rds-rb643-ukbiobank2/Data_Users/yh464/pheno/ukb')
_, ax = plt.subplots(2, 3, figsize = (30, 20))
ax = ax.reshape(6)
i = 0
for p in ['deg','degi','degc','clu','eff','mpl']:
    df = pd.read_table(f'{p}_local.txt', sep = '\s+')
    df = df.iloc[:,2:].describe()
    left = df.iloc[1:,:188].T
    right = df.iloc[1:, 188:].T
    left.columns = ['l_mean','l_std','l_min','l_q1','l_q2','l_q3','l_max']
    left.index = left.index.str.replace('L_','')
    right.columns = ['r_mean','r_std','r_min','r_q1','r_q2','r_q3','r_max']
    right.index = left.index
    df = pd.concat([left, right], axis = 1)
    
    ax[i].errorbar(x = df.l_q2, y = df.r_q2, xerr = abs(df[['l_q1','l_q3']].T-df.l_q2),
                   yerr = abs(df[['r_q1', 'r_q3']].T - df.r_q2), linestyle = '', markersize = 3)
    ax[i].set_title(p)
    i += 1

plt.savefig('diag.pdf', bbox_inches = 'tight')