# this script is for a single purpose and is not written using argparse
# change script to modify

fin = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/pheno/ukb/raw/ukb45629_psych.txt'
fout = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/pheno/ukb/ukb45629_psych'

import pandas as pd
import numpy as np
from fnmatch import fnmatch

dfin = pd.read_table(fin)
dfout = dfin.iloc[:,:2] # fid and iid
cols = dfin.columns
nsubj = dfin.shape[0]

# bipolar and MDD status
tmp = dfin['f.20126.0.0']
t = ((tmp==1) + (tmp==2)).astype(np.float16)
t[tmp.isna()] = np.nan
dfout['bip_online'] = t
t = (tmp>2).astype(np.float16)
t[tmp.isna()] = np.nan
dfout['mdd_online'] = t

# ever depressed for a whole week
tmp = [fnmatch(x, 'f.4598.*') for x in cols]
tmp = dfin.loc[:,tmp].values
t = np.full(nsubj,np.nan)
for i in range(nsubj):
  for j in range(tmp.shape[1]): # over instances
    if tmp[i,j] == 1: t[i] = 1 # true overwrites false
    if tmp[i,j] == 0: # false cannot overwrite true, so
      if t[i]**2 == 1: t[i] = -1
      else: t[i] = 0
t[t==-1] = np.nan
dfout['dep_whole_week'] = t

# antidepressants
tmp = [fnmatch(x, 'f.20546.*') for x in cols]
tmp = dfin.loc[:,tmp].values
t = np.full(nsubj,np.nan)
for i in range(nsubj):
  if any(~np.isnan(tmp[i,:])): t[i] = 0
  if any(tmp[i,:] == 1) or any(tmp[i,:] == 3): t[i] = 1
dfout['antidepressant'] = t

# ICD10 diagnoses
tmp = [fnmatch(x, 'f.41270.*') for x in cols]
tmp = dfin.loc[:,tmp]
t = np.full((nsubj,11),np.nan) # alz, mdd, scz, bip, adhd, ed, an, asd, ocd, anx, pd
for i in range(nsubj):
  if all(tmp.iloc[i,:].isna()): continue
  tr = tmp.iloc[i,:].values.astype('U')
  t[i,0] = (any([fnmatch(x, 'F00*') for x in tr]) or any([fnmatch(x, 'G30*') for x in tr]))
  t[i,1] = (any([fnmatch(x, 'F32*') for x in tr]) or any([fnmatch(x, 'F33*') for x in tr]))
  t[i,2] = any([fnmatch(x, 'F20*') for x in tr])
  t[i,3] = any([fnmatch(x, 'F31*') for x in tr])
  t[i,4] = any([fnmatch(x, 'F900') for x in tr])
  t[i,5] = any([fnmatch(x, 'F50*') for x in tr])
  t[i,6] = any([fnmatch(x, 'F500') for x in tr])
  t[i,7] = any([fnmatch(x, 'F84*') for x in tr])
  t[i,8] = any([fnmatch(x, 'F42*') for x in tr])
  t[i,9] = any([fnmatch(x, 'F41*') for x in tr])
  t[i,10]= any([fnmatch(x, 'F410') for x in tr])
dfout['adhd'] = t[:,4]
dfout['alz'] = t[:,0]
dfout['an'] = t[:,6]
dfout['anx'] = t[:,9]
dfout['bip'] = t[:,3]
dfout['ed'] = t[:,5]
dfout['mdd'] = t[:,1]
dfout['ocd'] = t[:,8]
dfout['pd'] = t[:,10]
dfout['scz'] = t[:,2]
dfout.to_csv(fout+'_binary.txt', sep = '\t', index=False)

# reset dfout for quantitated variables
dfout = dfin.iloc[:,:2] # fid and iid

# prospective memory
tmp = [fnmatch(x, 'f.20018.*') for x in cols]
tmp = dfin.loc[:,tmp].values
tmp[tmp==2] = .5 # score .5 points for correct recalls on 2nd attempt
dfout['prospective_memory'] = np.nanmean(tmp, axis = 1)

# educational attainment
tmp = [fnmatch(x, 'f.6138.*') for x in cols]
tmp = dfin.loc[:,tmp].values
tmp[tmp==-3] = np.nan # prefer not to answer
t = np.full(nsubj, np.nan)
for i in range(nsubj):
  tr = np.unique(tmp[i,:]) # union over all instances
  if 1 in tr: 
    t[i] = 16
  elif 6 in tr:
    t[i] = 14
  elif (5 in tr) or (2 in tr):
    t[i] = 12
  elif (3 in tr) or (4 in tr):
    t[i] = 9
  else: t[i] = 0
dfout['edu_years'] = t

# age completing education
tmp = [fnmatch(x, 'f.845.*') for x in cols] # 3 instances
tmp = dfin.loc[:,tmp].values
dfout['edu_age'] = np.nanmax(tmp, axis = 1) # max over instances

# fluid intelligence
tmp = [(fnmatch(x, 'f.20191.*') or fnmatch(x,'f.20016.*')) for x in cols]
tmp = dfin.loc[:,tmp].values
dfout['fluid'] = np.nanmean(tmp, axis = 1) # mean IQ over all instances

# neuroticism
dfout['neuroticism'] = dfin['f.20127.0.0']

# PHQ9
tmp = [(fnmatch(x, 'f.2051[4079183].*') or fnmatch(x,'f.2050[78].*')) for x in cols]
tmp = dfin.loc[:,tmp].values
tmp -= 1 # score using 0123 instead of 1234
t = tmp.sum(axis = 1)
t[t<0] = np.nan # "prefer not to answer" = -818
dfout['phq9'] = t

# GAD7
tmp = [(fnmatch(x, 'f.2051[256].*') or fnmatch(x,'f.2050[569].*') or 
        fnmatch(x, 'f.20520.*')) for x in cols]
tmp = dfin.loc[:,tmp].values
tmp -= 1 # score using 0123 instead of 1234
t = tmp.sum(axis = 1)
t[t<0] = np.nan # "prefer not to answer" = -818
dfout['gad7'] = t

dfout.to_csv(fout+'_quant.txt', sep = '\t', index = False)