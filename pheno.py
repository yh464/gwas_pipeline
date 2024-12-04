#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1.0: 2023-07-02
Version 1.1: 2024-11-26

Generates graph metrics for connectomes

Preceding workflow:
    fmri_parcel_batch.py
Requires following inputs: 
    connectivity matrices after wavelets filtering

Changelog:
    changed output format so that index = phenotype name, columns = phenotype group name
'''

import numpy as np
import scipy.stats as sts
import pandas as pd
import bct
import argparse as ap
import os
import time

# input argument processing
parser = ap.ArgumentParser(description='This programme processes the connectome '+
                           ' for one single individual for imaging derived phenotypes')
parser.add_argument('subj',help = 'Subject ID')
parser.add_argument('-i','--in',dest = '_in', help =
    'Target file to screen',
    default = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Imaging/'+
    '%sub/func/fMRI/parcellations/HCP.fsaverage.aparc_seq/Connectivity_sc2345.txt')
parser.add_argument('-o','--out',dest = 'out', help = 'Output directory',
                    default = '../pheno/ukb/')
parser.add_argument('-f','--force', dest = 'force', help = 'Force output',
                    default = False,action = 'store_true')
args = parser.parse_args()
args.out = os.path.realpath(args.out)

# nroi: HCP = 376, 500sym = 334, aparc = 84, economo = 102, sjh = 1027
nroi = 376

# specify directories
in_filename = args._in.replace('%sub',args.subj)

if not os.path.isdir(args.out): os.system(f'mkdir -p {args.out}/')              # creates parent folders
if not os.path.isdir(f'{args.out}/global/'): os.mkdir(f'{args.out}/global/')
if not os.path.isdir(f'{args.out}/local/'): os.mkdir(f'{args.out}/local/')
if not os.path.isdir(f'{args.out}/global_asym/'): os.mkdir(f'{args.out}/global_asym/')
if not os.path.isdir(f'{args.out}/local_asym/'): os.mkdir(f'{args.out}/local_asym/')

if not os.path.isfile(in_filename):
  raise ValueError('No connectome found for the subject')

# loading the connectome  
tic = time.perf_counter()
rsc = np.loadtxt(in_filename,delimiter = ',')
if nroi != rsc.shape[0]:
  raise ValueError('Connectome matrix does not have the right shape')
order = np.concatenate((np.arange(0,8),np.arange(16,196),
                        np.arange(8,16),np.arange(196,376)))
rsc = rsc[order,:][:,order]                                                    # re-order to 188 L, 188 R   
nroi_5 = int(nroi/2)                                                           # force int for indexing
if np.isnan(rsc).any():
  raise ValueError('NaN connectivity')
  
# node labels
nodes = np.loadtxt(
  '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/params/HCP.fsaverage_4mm_names.txt', 
  dtype ='U')
nodes_u = nodes[:int(nodes.size/2)].copy()
for i in range(nodes_u.size):
  nodes_u[i] = nodes_u[i][2:] # unilateral label

# check progress
skip = True
if args.force: skip = False

if skip:
  try:
    tmp = np.loadtxt(f'{args.out}/global/{args.subj}.txt')
    if tmp.size != 17: skip = False
  except: skip = False

if skip:
  try:
    tmp = np.loadtxt(f'{args.out}/global_asym/{args.subj}.txt')
    if tmp.size != 17: skip = False
  except: skip = False

if skip:
  try:
    tmp = pd.read_csv(f'{args.out}/local/{args.subj}.txt')
    if tmp.shape[0] != nroi or tmp.shape[1] != 7: skip = False
  except: skip = False

if skip:
  try:
    tmp = pd.read_csv(f'{args.out}/local_asym/{args.subj}.txt')
    if tmp.shape[0] != nroi_5 or tmp.shape[1] != 21: skip = False
  except: skip = False

if skip:
  raise ValueError('Stats already generated')

# set negative correlations and autocorrelations to zero
rsc[rsc<0] = 0
for i in range(nroi): rsc[i,i] = 0

# Fisher transform
rsc = np.arctanh(rsc)
rsc = rsc/2 + rsc.T/2                                                          # enforce symmetry
rsc = rsc/rsc.max()                                                            # re-scaling as recommended by BCT

# neighbourhood
rsc_bin = rsc>0
deg_bin = rsc_bin.sum(axis=1)

# degree measures
deg_local = rsc.mean(axis = 1)                                                 # out by nroi on purpose
deg_global = deg_local.mean()
degi_local = np.concatenate([rsc[:nroi_5,:nroi_5].mean(axis=1),
              rsc[nroi_5:,nroi_5:].mean(axis=1)])
degi_global = degi_local.mean()
degi_hem = np.array([degi_local[:nroi_5].mean(),degi_local[nroi_5:].mean()])
degc_local = np.concatenate([rsc[:nroi_5,nroi_5:].mean(axis=1),
              rsc[nroi_5:,:nroi_5].mean(axis=1)])
degc_global = degc_local.mean()

# path length measures - needed for later use in random models
def efficiency(rsc):
  rsc_len = rsc.copy()
  rsc_len[rsc_bin] = 1/rsc_len[rsc_bin]
  dijkstra,_ = bct.distance_wei(rsc_len)
  dijkstra_inv = dijkstra.copy()                                               # invert Dijkstra distance
  dijkstra_inv[dijkstra_inv > 0] = 1/dijkstra_inv[dijkstra_inv > 0]
  mpl_local = dijkstra.mean(axis=1)
  mpl_global = mpl_local.mean()
  mpl_hem = np.array([dijkstra[:nroi_5,:nroi_5].mean(),                        # intrahemispheric
             dijkstra[nroi_5:,nroi_5:].mean()])
  
  # efficiency measures
  dijkstra_3,_ = bct.distance_wei(rsc_len**(1/3))
  dijkstra_3[dijkstra_3>0] = 1/dijkstra_3[dijkstra_3>0]                        # as per bct
  eff_global = dijkstra_inv.sum()/nroi/(nroi-1)                                # there is a need to decompose the code to reduce repeated computation! 
  eff_hem = np.array([dijkstra_inv[:nroi_5,:nroi_5].sum()/nroi_5/(nroi_5-1),   # intrahemispheric efficiency 
             dijkstra_inv[nroi_5:,nroi_5:].sum()/nroi_5/(nroi_5-1),])
  eff_local = np.zeros(nroi)
  for i in range(nroi):
    n_i = rsc_bin[:,i]
    #num_i = (rsc[n_i,:][:,n_i]**(1/3)@rsc[n_i,:][:,n_i]**(1/3))
    num_i = (rsc[n_i,i]**(1/3)).reshape((deg_bin[i],1))
    num_i = (num_i@num_i.T) * dijkstra_3[n_i,:][:,n_i]                         # NB bctpy uses an older version
    eff_local[i] = num_i.sum()/deg_bin[i]/(deg_bin[i]-1)
  return mpl_global, mpl_hem, mpl_local, eff_global, eff_hem, eff_local

mpl_global, mpl_hem, mpl_local, eff_global, eff_hem, eff_local = efficiency(rsc)

# clustering measures
clu_local = bct.clustering_coef_wu(rsc)
clu_global = clu_local.mean()
clu_hem = np.array([bct.clustering_coef_wu(rsc[:nroi_5,:nroi_5]).mean(),       # intrahemispheric
           bct.clustering_coef_wu(rsc[nroi_5:,nroi_5:]).mean()])

# small_worldness measures
smw_global = clu_global/mpl_global
smw_hem = clu_hem/mpl_hem

# centrality measures
bet_local = bct.betweenness_wei(rsc)/(nroi-1)/(nroi-2)
toc = time.perf_counter()-tic
print(f'finished global and regional phenotypes. time = {toc:.3f} seconds')

# hemispheric asymmetry measures - intrahemispheric networks
degi_asym_abs = degi_hem[0]-degi_hem[1]
degi_asym_frac = 2*degi_asym_abs/degi_hem.sum()                                # difference / mean

mpl_asym_abs = mpl_hem[0] - mpl_hem[1]
mpl_asym_frac = 2*mpl_asym_abs/mpl_hem.sum()

eff_asym_abs = eff_hem[0] - eff_hem[1]
eff_asym_frac = 2*eff_asym_abs/eff_hem.sum()

clu_asym_abs = clu_hem[0] - clu_hem[1]
clu_asym_frac = 2*clu_asym_abs/clu_hem.sum()

smw_asym_abs = smw_hem[0] - smw_hem[1]
smw_asym_frac = 2*smw_asym_abs/smw_hem.sum()

# hemispheric asymmetry measures - correlational
def squareform(mat):
    tmp = []
    for i in range(mat.shape[0]):
        tmp.append(mat[i,:i])
    tmp = np.concatenate(tmp)
    return tmp

deg_asym_corr = sts.spearmanr(deg_local[:nroi_5],deg_local[nroi_5:])           # this is in essence a symmetry measure
deg_asym_corr = -np.arctanh(deg_asym_corr[0])                                     # normalise to z score, -ve to reflect asymmetry

degi_asym_corr = sts.spearmanr(degi_local[:nroi_5],degi_local[nroi_5:])          
degi_asym_corr = -np.arctanh(degi_asym_corr[0])

degc_asym_corr = sts.spearmanr(degc_local[:nroi_5],degc_local[nroi_5:])          
degc_asym_corr = -np.arctanh(degc_asym_corr[0])

mpl_asym_corr = sts.spearmanr(mpl_local[:nroi_5],mpl_local[nroi_5:])          
mpl_asym_corr = -np.arctanh(mpl_asym_corr[0])

eff_asym_corr = sts.spearmanr(eff_local[:nroi_5],eff_local[nroi_5:])          
eff_asym_corr = -np.arctanh(eff_asym_corr[0])

clu_asym_corr = sts.spearmanr(clu_local[:nroi_5],clu_local[nroi_5:])          
clu_asym_corr = -np.arctanh(clu_asym_corr[0])

bet_asym_corr = sts.spearmanr(bet_local[:nroi_5],bet_local[nroi_5:])          
bet_asym_corr = -np.arctanh(bet_asym_corr[0])

# asymmetry of all intra-hemispheric connectivities
con_asym_corr = sts.spearmanr(squareform(rsc[:nroi_5,:nroi_5]), squareform(rsc[nroi_5:, nroi_5:]))
con_asym_corr = -np.arctanh(con_asym_corr[0])

# local asymmetry measures
deg_nasym_abs = deg_local[:nroi_5]-deg_local[nroi_5:]
deg_nasym_frac = 2*deg_nasym_abs/(deg_local[:nroi_5]+deg_local[nroi_5:])
deg_nasym_frac[np.isnan(deg_nasym_frac)]=0  
deg_nasym_rank = sts.rankdata(deg_local)                                       # ties = average rank
deg_nasym_rank = (deg_nasym_rank[:nroi_5]-deg_nasym_rank[nroi_5:])/nroi        # normalise rank difference to (0,1)

degi_nasym_abs = degi_local[:nroi_5]-degi_local[nroi_5:]
degi_nasym_frac = 2*degi_nasym_abs/(degi_local[:nroi_5]+degi_local[nroi_5:])
degi_nasym_frac[np.isnan(degi_nasym_frac)]=0  
degi_nasym_rank = sts.rankdata(degi_local)
degi_nasym_rank = (degi_nasym_rank[:nroi_5]-degi_nasym_rank[nroi_5:])/nroi

degc_nasym_abs = degc_local[:nroi_5]-degc_local[nroi_5:]
degc_nasym_frac = 2*degc_nasym_abs/(degc_local[:nroi_5]+degc_local[nroi_5:])
degc_nasym_frac[np.isnan(degc_nasym_frac)]=0  
degc_nasym_rank = sts.rankdata(degc_local)
degc_nasym_rank = (degc_nasym_rank[:nroi_5]-degc_nasym_rank[nroi_5:])/nroi

mpl_nasym_abs = mpl_local[:nroi_5]-mpl_local[nroi_5:]
mpl_nasym_frac = 2*mpl_nasym_abs/(mpl_local[:nroi_5]+mpl_local[nroi_5:])
mpl_nasym_frac[np.isnan(mpl_nasym_frac)]=0  
mpl_nasym_rank = sts.rankdata(mpl_local)
mpl_nasym_rank = (mpl_nasym_rank[:nroi_5]-mpl_nasym_rank[nroi_5:])/nroi

eff_nasym_abs = eff_local[:nroi_5]-eff_local[nroi_5:]
eff_nasym_frac = 2*eff_nasym_abs/(eff_local[:nroi_5]+eff_local[nroi_5:])
eff_nasym_frac[np.isnan(eff_nasym_frac)]=0  
eff_nasym_rank = sts.rankdata(eff_local)
eff_nasym_rank = (eff_nasym_rank[:nroi_5]-eff_nasym_rank[nroi_5:])/nroi

clu_nasym_abs = clu_local[:nroi_5]-clu_local[nroi_5:]
clu_nasym_frac = 2*clu_nasym_abs/(clu_local[:nroi_5]+clu_local[nroi_5:])
clu_nasym_frac[np.isnan(clu_nasym_frac)]=0  
clu_nasym_rank = sts.rankdata(clu_local)
clu_nasym_rank = (clu_nasym_rank[:nroi_5]-clu_nasym_rank[nroi_5:])/nroi

bet_nasym_abs = bet_local[:nroi_5]-bet_local[nroi_5:]
bet_nasym_frac = 2*bet_nasym_abs/(bet_local[:nroi_5]+bet_local[nroi_5:])
bet_nasym_frac[np.isnan(bet_nasym_frac)]=0                                    
bet_nasym_rank = sts.rankdata(bet_local)
bet_nasym_rank = (bet_nasym_rank[:nroi_5]-bet_nasym_rank[nroi_5:])/nroi

toc = time.perf_counter()-tic
print(f'finished asymmetry phenotypes. time = {toc:.3f} seconds')

# output files
# global measures
out_filename = f'{args.out}/global/{args.subj}.txt'
# output order: 7 global measures and 5 hemispheric measures in blocks of 2
df = pd.DataFrame(dict(
    deg_global = deg_global, degi_global = degi_global, degc_global = degc_global,
    mpl_global = mpl_global, eff_global = eff_global, clu_global = clu_global,
    smw_global = smw_global,
    degi_l = degi_hem[0], degi_r = degi_hem[1], mpl_l = mpl_hem[0], mpl_r = mpl_hem[1],
    eff_l = eff_hem[0], eff_r = eff_hem[1], clu_l = clu_hem[0], clu_r = clu_hem[1],
    smw_l = smw_hem[0], smw_r = [smw_hem[1]]
    ), index = ['global']).T
df.to_csv(out_filename, sep ='\t', index = True, header = True)

# local measures      
out_filename = f'{args.out}/local/{args.subj}.txt'
df = pd.DataFrame(dict(deg_local = deg_local, degc_local = degc_local, 
  degi_local = degi_local, mpl_local = mpl_local, eff_local = eff_local,
  clu_local = clu_local, bet_local = bet_local), index = nodes)
df.to_csv(out_filename,sep='\t', index = True, header = True)

# global asymmetries
out_filename = f'{args.out}/global_asym/{args.subj}.txt'
df = pd.DataFrame(dict(degi_asym_abs = degi_asym_abs, degi_asym_frac = degi_asym_frac,
                       mpl_asym_abs = mpl_asym_abs, mpl_asym_frac = mpl_asym_frac,
                       eff_asym_abs = eff_asym_abs, eff_asym_frac = eff_asym_frac,
                       clu_asym_abs = clu_asym_abs, clu_asym_frac = clu_asym_frac,
                       smw_asym_abs = smw_asym_abs, smw_asym_frac = smw_asym_frac,
                       deg_asym_corr = deg_asym_corr, degi_asym_corr = degi_asym_corr,
                       degc_asym_corr = degc_asym_corr, mpl_asym_corr = mpl_asym_corr,
                       eff_asym_corr = eff_asym_corr, clu_asym_corr = clu_asym_corr,
                       bet_asym_corr = bet_asym_corr, con_asym_corr = [con_asym_corr]
                       ), index = ['global_asym']).T
df.to_csv(out_filename, sep = '\t', index = True, header = True)

# local asymmetries
out_filename = f'{args.out}/local_asym/{args.subj}.txt'
df = pd.DataFrame(dict(
  deg_nasym_abs = deg_nasym_abs, deg_nasym_frac = deg_nasym_frac, deg_nasym_rank = deg_nasym_rank,
  degi_nasym_abs = degi_nasym_abs, degi_nasym_frac = degi_nasym_frac, degi_nasym_rank = degi_nasym_rank,
  degc_nasym_abs = degc_nasym_abs, degc_nasym_frac = degc_nasym_frac, degc_nasym_rank = degc_nasym_rank,
  mpl_nasym_abs = mpl_nasym_abs, mpl_nasym_frac = mpl_nasym_frac, mpl_nasym_rank = mpl_nasym_rank,
  eff_nasym_abs = eff_nasym_abs, eff_nasym_frac = eff_nasym_frac, eff_nasym_rank = eff_nasym_rank,
  clu_nasym_abs = clu_nasym_abs, clu_nasym_frac = clu_nasym_frac, clu_nasym_rank = clu_nasym_rank,
  bet_nasym_abs = bet_nasym_abs, bet_nasym_frac = bet_nasym_frac, bet_nasym_rank = bet_nasym_rank), index = nodes_u)
df.to_csv(out_filename,sep='\t', index = True, header = True)