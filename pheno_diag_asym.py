#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2024-11-27

This is a diagnostic script to understand asymmetry phenotypes of the network
Intended to be one-off
'''

import argparse
parser = argparse.ArgumentParser(description = 'FDR script for local spin perm tests')
parser.add_argument('-f', dest = 'force', action = 'store_true', default = False,
                    help = 'force overwrite')
args = parser.parse_args()
force = args.force

import numpy as np
import pandas as pd
import statsmodels.api as sm
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
from _plots import corr_heatmap

# aesthetics
cdict = dict(red = ((0,0,0),(1/2,1,1),(1,.8,.8)),
             green = ((0,0,0),(1/2,1,1),(1,0,0)),
             blue = ((0,.8,.8),(1/2,1,1),(1,0,0)))
cmap_name = 'redblue'
cmap = mpl.colors.LinearSegmentedColormap(cmap_name,cdict,1024)
try:
  mpl.colormaps.register(cmap)
except:
  pass
sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 2})
sns.set_theme(style = 'ticks')

os.chdir('/rds/user/yh464/rds-rb643-ukbiobank2/Data_Users/yh464/pheno/ukb')
if not os.path.isdir('diagnostics'): os.mkdir('diagnostics')

# distribution of global phenotypes
glob = pd.read_table('global_graph.txt', sep = '\\s+').set_index('IID')
asym = pd.read_table('global_asym_diag.txt', sep = '\\s+').set_index('IID')
rois = pd.Index(open('deg_local.txt').readline().replace('\n','').split()[2:]) # only read header
if not os.path.isfile('diagnostics/global_asym_dist.svg') or force:
    tmp = glob[['deg_global','degi_global','degc_global','eff_global','mpl_global',
                'clu_global','smw_global']]
    tmp.columns = ['Degree','Ipsilateral\ndegree','Contralateral\ndegree','Efficiency',
                   'Path length', 'Clustering', 'Small-worldness']
    tmp = tmp.melt(var_name = 'phenotype')
    plt.subplots(figsize = (5,5))
    sns.histplot(tmp, x = 'value', hue = 'phenotype', element = 'step')
    plt.title('Global phenotypes', fontsize = 'large')
    plt.savefig('diagnostics/global_graph_dist.pdf', bbox_inches = 'tight')
    plt.close()
    tmp = asym[['deg_asym_corr','degi_asym_corr','degc_asym_corr','eff_asym_corr',
                'mpl_asym_corr','clu_asym_corr']]
    tmp.columns = ['Degree','Ipsilateral\ndegree','Contralateral\ndegree','Efficiency',
                   'Path length', 'Clustering']
    tmp = tmp.melt(var_name = 'phenotype')
    plt.subplots(figsize = (5,5))
    sns.histplot(tmp, x = 'value', hue = 'phenotype', element = 'step')
    plt.title('Asymmetry phenotypes', fontsize = 'large')
    plt.savefig('diagnostics/global_asym_dist.pdf', bbox_inches = 'tight')
    plt.close()

# correlate left v right
fig0, ax0 = plt.subplots(2, 3, figsize = (15, 10))
ax0 = ax0.reshape(6)
fig1, ax1 = plt.subplots(2, 3, figsize = (15, 10))
ax1 = ax1.reshape(6)
i = 0
for p, heading in zip(['deg','degi','degc','clu','eff','mpl'],
                      ['Degree','Ipsilateral degree','Contralateral degree',
                       'Efficiency','Path length', 'Clustering']):
    tmp = asym[f'{p}_asym_corr']
    asym_max = asym.loc[tmp.idxmax(), 'FID']
    asym_min = asym.loc[tmp.idxmin(), 'FID']
    asym_q1 = asym.loc[tmp==tmp.quantile(.25, 'nearest'),'FID']
    asym_q2 = asym.loc[tmp==tmp.quantile(.5, 'nearest'),'FID']
    asym_q3 = asym.loc[tmp==tmp.quantile(.75, 'nearest'),'FID']
    
    df = pd.read_table(f'{p}_local.txt', sep = '\\s+').set_index('FID').iloc[:,1:]
    # for x, lab in zip([asym_min, asym_q1, asym_q2, asym_q3, asym_max],
    #                   ['min','25%','median','75%','max']):
    # for x, lab in zip([asym_q1, asym_q2, asym_q3],
    #                   ['25%','median','75%']):
    #     tmp = df.loc[x,:]
    #     if tmp.size > 376:
    #         tmp = tmp.iloc[0,:]
    #     tmp = tmp.to_numpy().reshape(376)
    #     ax0[i].scatter(tmp[:188], tmp[188:], s=2, label = f'{lab} asymmetry')
    ax0[i].scatter(df.iloc[21894,:188],df.iloc[21894,188:], s = 3, label = 'Example individual')
    df = df.describe()
    left = df.iloc[1:,:188].T
    right = df.iloc[1:, 188:].T
    left.columns = ['l_mean','l_std','l_min','l_q1','l_q2','l_q3','l_max']
    left.index = left.index.str.replace('L_','')
    right.columns = ['r_mean','r_std','r_min','r_q1','r_q2','r_q3','r_max']
    right.index = left.index
    df = pd.concat([left, right], axis = 1)
    
    if not os.path.isfile('diagnostics/region_LvsR.pdf') or force:
        ax0[i].errorbar(x = df.l_q2, y = df.r_q2, xerr = abs(df[['l_q1','l_q3']].T-df.l_q2),
                       yerr = abs(df[['r_q1', 'r_q3']].T - df.r_q2), linestyle = '', markersize = 2,
                       label = 'Quartiles', color = 'k', linewidth = 0.4,zorder = 0.5)
        ax0[i].set_title(heading, fontsize = 'large')
        ax0[i].axline((0,0), slope = 1, color = 'k')
        # ax0[i].legend(fontsize = 'large')
        ax0[i].set_xlabel('L hemisphere', fontsize = 'large')
        ax0[i].set_ylabel('R hemisphere', fontsize = 'large')
    
    if not os.path.isfile('diagnostics/asym_compare.png') or force:
        try: ax1[i].scatter(x = asym[f'{p}_asym_corr'], y = asym[f'{p}_asym_pear'], label = 'pear')
        except: pass
        try: ax1[i].scatter(x = asym[f'{p}_asym_corr'], y = asym[f'{p}_asym_frac'], label = 'frac')
        except: pass
        try: ax1[i].scatter(x = asym[f'{p}_asym_corr'], y = asym[f'{p}_asym_abs'], label = 'abs')
        except: pass
        ax1[i].set_title(heading)
        ax1[i].legend()
    i += 1

if not os.path.isfile('diagnostics/region_LvsR.pdf') or force:
    fig0.savefig('diagnostics/region_LvsR.pdf', bbox_inches = 'tight')
if not os.path.isfile('diagnostics/asym_compare.png') or force:
    fig1.savefig('diagnostics/asym_compare.png', bbox_inches = 'tight')
plt.close()

# randomly take 5000 individuals
rng = np.random.Generator(np.random.MT19937(114514))
eid = rng.choice(glob['FID'], size = 5000, replace = False)
eid = sorted(eid)


for roi in ['V3A','PF','Ig','MT','pOFC']:
  if not os.path.isfile(f'diagnostics/regional_v_asym_{roi}.pdf') or force:
    # random region L/R difference v asymmetry
    df = pd.read_table('deg_local.txt', sep = '\\s+').set_index('FID').iloc[:,1:]
    _, ax = plt.subplots(1,3,figsize = (15,5))
    df1 = df.loc[eid,:]
    asym_test = asym.loc[asym.FID.isin(eid),['FID','deg_asym_corr']]
    df1 = pd.merge(df1, asym_test, on = 'FID').set_index('FID').sort_values(
        by = f'L_{roi}_ROI')
    sns.scatterplot(df1, x = f'L_{roi}_ROI', y = f'R_{roi}_ROI', hue = 'deg_asym_corr',
                    s = 1, linewidths = 0, palette = 'redblue', ax = ax[0])
    
    df1['L - R'] = df1[f'L_{roi}_ROI']- df1[f'R_{roi}_ROI']
    sns.scatterplot(df1, x = 'deg_asym_corr', y = 'L - R', hue = f'L_{roi}_ROI',
        s = 2, linewidths = 0, palette = 'redblue', ax = ax[1])
    
    df1 = df1.melt(id_vars = 'deg_asym_corr', value_vars = [f'L_{roi}_ROI',f'R_{roi}_ROI'], 
                   var_name = 'hemisphere')
    sns.scatterplot(df1, x = 'deg_asym_corr', y = 'value', hue = 'hemisphere',
        s = 2, linewidths = 0, ax = ax[2])
    # ax[2].set_title('Example region: {roi}')
    plt.savefig(f'diagnostics/regional_v_asym_{roi}.pdf', bbox_inches = 'tight')
    plt.close()

# split-half regression
top_half = pd.DataFrame(columns = ['deg','degi','degc','clu','eff','mpl'],
                     index = rois, data = np.nan)
bottom_half = top_half.copy()
betas = top_half.copy()
diff_corr = top_half.copy()
for p in ['deg','degi','degc','clu','eff','mpl']:
    df = pd.read_table(f'{p}_local.txt', sep = '\\s+').set_index('FID').iloc[:,1:]
    lrtotal = pd.DataFrame(columns = rois[:188].str.replace('L_','') + '_total',
        index = df.index, data = df.iloc[:,:188].to_numpy() + df.iloc[:,188:376].to_numpy())
    lrdiff = pd.read_table(f'{p}_nasym_abs.txt', sep='\\s+').set_index('FID').iloc[:, 1:]
    lrdiff.columns = rois[:188].str.replace('L_','') + '_diff'
    df = pd.concat([df, lrtotal, lrdiff], axis = 1)
    tmp = pd.concat([asym[f'{p}_asym_corr'], glob[f'{p}_global']], axis = 1)
    df = pd.concat([df, tmp], axis = 1)
    for roi in rois:
        total = roi.replace('L_','').replace('R_','') + '_total'
        diff = roi.replace('L_','').replace('R_','') + '_diff'
        m = np.median(df[total])
        top_half.loc[roi, p] = df.loc[df[total]>m,roi].corr(
            df.loc[df[total]>m,f'{p}_asym_corr'])
        bottom_half.loc[roi, p] = df.loc[df[total]<m,roi].corr(
            df.loc[df[total]<m,f'{p}_asym_corr'])
        diff_corr.loc[roi,p] = df[f'{p}_asym_corr'].corr(df[diff])
        try:
            tmp = df[[f'{p}_asym_corr',total, roi]].dropna()
            tmp['contra'] = tmp[total]-tmp[roi]
            md = sm.OLS(endog = tmp[f'{p}_asym_corr'],
                    # exog = tmp[['contra', roi]]).fit()
                    exog = tmp[[total, roi]]).fit()
            betas.loc[roi, p] = md.params[roi]
        except: pass
        # if md.pvalues[roi] < 0.05: betas_sigonly.loc[roi,p] = md.params[roi]

betas.to_csv('diagnostics/adjpcorr_adjL+R_regional_asym.txt', sep ='\t', index = True)
top_half.to_csv('diagnostics/adjpcorr_tophalf_L+R_regional_asym.txt', sep ='\t', index = True)
bottom_half.to_csv('diagnostics/adjpcorr_bottomhalf_L+R_regional_asym.txt', sep ='\t', index = True)
diff_corr.to_csv('diagnostics/pcorr_nasym_asym.txt', sep = '\t', index = True)
# # correlate local with asymmetry and global
# glob_corr = []
# asym_corr = []
# pear_corr = []
# frac_corr = []
# for p in ['deg','degi','degc','clu','eff','mpl']:
#     df = pd.read_table(f'{p}_local.txt', sep = '\\s+')
#     rois = df.columns[2:] # FID IID ...
#     glob_temp = glob[['IID',f'{p}_global']]
#     df = pd.merge(glob_temp, df)
#     if not p in ['deg','degc']:
#         asym_temp = asym[['IID',f'{p}_asym_corr', f'{p}_asym_pear', f'{p}_asym_frac']]
#     else:
#         asym_temp = asym[['IID',f'{p}_asym_corr', f'{p}_asym_pear']]
#     df = pd.merge(asym_temp, df)
    
#     cor = df.drop(['FID','IID'], axis = 1).corr()
#     glob_corr.append(cor.loc[rois, f'{p}_global'])
#     asym_corr.append(cor.loc[rois, f'{p}_asym_corr'])
#     pear_corr.append(cor.loc[rois, f'{p}_asym_pear'])
#     if not p in ['deg','degc']:
#         frac_corr.append(cor.loc[rois, f'{p}_asym_frac'])

# glob_corr = pd.concat(glob_corr, axis = 1)
# glob_corr.index.name = 'label'
# glob_corr.to_csv('diagnostics/pcorr_regional_global.txt', sep = '\t', index = True)
# asym_corr = pd.concat(asym_corr, axis = 1)
# asym_corr.index.name = 'label'
# asym_corr.to_csv('diagnostics/pcorr_regional_asym_spearman.txt', sep = '\t', index = True)
# pear_corr = pd.concat(pear_corr, axis = 1)
# pear_corr.index.name = 'label'
# pear_corr.to_csv('diagnostics/pcorr_regional_asym_pearson.txt', sep = '\t', index = True)
# frac_corr = pd.concat(frac_corr, axis = 1)
# frac_corr.index.name = 'label'
# frac_corr.to_csv('diagnostics/pcorr_regional_asym_frac.txt', sep = '\t', index = True)