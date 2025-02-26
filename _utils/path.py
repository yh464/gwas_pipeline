# -*- coding: utf-8 -*-
#!/usr/env/bin python3

'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2025-01-20

This utility manages the nomenclature of file names / paths
'''

import os
import numpy as np
import re
import pandas as pd

class normaliser():
    def __init__(self, _dir = os.path.realpath('../path/'), _dict = 'dict.txt'):
        
        self._dict_file = f'{_dir}/{_dict}'
        
        # initialise file
        if not os.path.isdir(_dir): os.mkdir(_dir)
        if not os.path.isfile(self._dict_file):
            df = pd.DataFrame(dict(
                before = ['pfdr','fdr','pheno','stderr','logor'],
                after = ['q','q','phenotype','se','beta']
                )).set_index('before')
            df.to_csv(self._dict_file, sep = '\t')
        
        # load file
        self._dict = pd.read_table(self._dict_file, index_col = 'before').fillna('')
    
    def load(self):
        self._dict = pd.read_table(self._dict_file, index_col = 'before').fillna('')
    
    def save(self):
        self._dict.to_csv(self._dict_file, sep = '\t', index = True, header = True)
    
    def append(self, before, after):
        # accepts new entries only, does not overwrite existing entries
        try:
            _ = self._dict.loc[before,'after']
            print('WARNING: cannot overwrite existing dictionary entry')
        except:
            self._dict = pd.concat((self._dict, 
                pd.DataFrame(index = [before], columns =['after'], data = after)))
            self.save()
    
    def update(self, before, after):
        # updates existing entries and appends new entries
        try:
            self._dict.loc[before,'after'] = after
            self.save()
            print(f'WARNING: overwriting existing entry: {before}')
        except:
            self._dict = pd.concat((self._dict, 
                pd.DataFrame(index = [before], columns =['after'], data = after)))
            self.save()
    
    def _normalise_df(self, df_in):
        def n(series, x, y):
            series = pd.Series(series)
            if series.dtype in [int, float]: return series
            if x[0] == '_':
                series = series.str.replace(x, y, regex = True, case = False)
            else:
                series = series.replace(x,y).replace(x.lower(), y).replace(x.upper(),y)
                series = '_' + series + '_'
                series = series.str.replace(f'_{x}_',f'_{y}_', regex = True, case = False).str.replace(
                f' {x}_',f' {y}_', regex = True, case = False).str.removeprefix('_').str.removesuffix('_')
            return series
        
        df = df_in.copy()
        for x in self._dict.index:
            y = str(self._dict.loc[x, 'after'])
            for col in df.columns:
                if col == 'SNP': continue # NEVER change SNP IDs
                try: df.loc[:,col] = n(df[col], x, y)
                except: pass
                
            col = df.columns.to_frame()
            for c in col.columns:
                try: col.loc[:,c] = n(col[c],x,y)
                except: pass
            if col.shape[1] > 1:
                df.columns = pd.MultiIndex.from_frame(col)
            else: df.columns = col.iloc[:,0]
            idx = df.index.to_frame()
            for c in idx.columns:
                try: idx.loc[:,c] = n(idx[c],x,y)
                except: pass
            if idx.shape[1] > 1:
                df.set_index(pd.MultiIndex.from_frame(idx), inplace = True)
            else: df.index = idx.iloc[:,0]
        return df
    
    def normalise(self, data, backup = None):
        # read table if data is a file
        if type(data) == str and os.path.isfile(data):
            if backup == None: backup = f'{data}.bak'
            os.system(f'cp {data} {backup}')
            if data[-3:] == 'csv':                
                df = pd.read_csv(data)
                df = self._normalise_df(df)
                df.to_csv(data, index = False)
            else: 
                df = pd.read_table(data, sep = '\\s+')
                df = self._normalise_df(df)
                df.to_csv(data, index = False, sep = '\t')
            return df
        
        # process input DataFrame object
        elif type(data) == pd.DataFrame:
            df = data
        # or try to coerse input object into pd.DataFrame
        else: df = pd.DataFrame(data)
        
        return self._normalise_df(df)

class project():
    def __init__(self,
                 _dir = '../path/', # uses relative path, so defaults to the relative path to the wd
                 _dict = 'wildcards.txt', # placeholders of a certain format
                 _flow = 'workflow.txt', # workflow, shows file name, input and output
                 ):
        
        # file names
        self._dict_file = _dir+_dict
        self._flow_file = _dir+_flow
        
        # create files with header
        if not os.path.isdir(_dir):
            os.mkdir(_dir)
        if not os.path.isfile(self._dict_file):
            f = open(self._dict_file,'w')
            f.write('var\tformat\tdescription')
            f.close()
            del f
        if not os.path.isfile(self._flow_file):
            f = open(self._flow_file,'w')
            f.write('path\toutput from\tinput to')
            f.close()
            del f
        
        # load files
        self._dict = np.loadtxt(self._dict_file, dtype = '<U1024', delimiter = '\t')
        if len(self._dict.shape) == 1:
            self._dict = self._dict.reshape((1,self._dict.size))
        self._flow = np.loadtxt(self._flow_file, dtype = '<U1024', delimiter = '\t')
        if len(self._flow.shape) == 1:
            self._flow = self._flow.reshape((1,self._flow.size))
    
    def add_var(self,name, fmt, desc):
        # format variable name
        name = str(name)
        if ord(name[0]) != ord('%'):
            name = '%' + name
        
        # sanity check
        for i in self._dict:
            if name == i[0]:
                if fmt == i[1] and desc == i[2]:
                    return
                else:
                    raise ValueError('Var name already occupied')
        # fmt must be a valid Regular Expression /lib/re
        
        # save file
        self._dict = np.vstack((self._dict, [name, fmt, desc]))
        np.savetxt(self._dict_file, self._dict, delimiter = '\t', fmt = '%s')
    
    def sanity_check(self, file_name, template_path):
        # unfortunately it is not possible to re-create the template path from file names
        # because multiple variables may be of the same format
        # template path contains the above variables like %subj
        for i in self._dict[1:,:]:# iterates over an entire row
            template_path = template_path.replace(i[0],i[1]) # name and fmt of _dict
        res = re.search(template_path, file_name)
        if type(res) == type(None):
            return False # may also raise an error for incorrect file names and return None for correct names
        else: return True
    
    def add_input(self, file, script):
        # fail-safe
        file = os.path.relpath(file) # uses relative path to the wd, i.e. ***/scripts/
        script = os.path.relpath(script)
        files = self._flow[:,0]
        
        # if the file is already included in the workflow file
        if len(np.argwhere(files==file)) == 1:
            idx = np.argwhere(files==file)[0][0]
            scripts = self._flow[idx,-1].split(', ')
            if len(scripts) == 0: self._flow[idx,-1] = script
            elif script in scripts: pass # do nothing if this input is already recorded
            else: self._flow[idx,-1] += f', {script}'
        # if the file is not otherwise found in the workflow file
        elif len(np.argwhere(files==file)) == 0:
            self._flow = np.vstack((self._flow,
                np.array([file,'',script], dtype = '<U1024')))
        else: raise ValueError('Check workflow file, repetitive entries found')
        
        # save file
        np.savetxt(self._flow_file, self._flow, delimiter = '\t', fmt = '%s')
    
    def add_output(self, file, script): # script means the script that generates the file
        # fail-safe
        file = os.path.relpath(file) # uses relative path to the wd, i.e. ***/scripts/
        script = os.path.relpath(script)
        files = self._flow[:,0]
        
        # if the file is already included in the workflow file
        if len(np.argwhere(files==file)) == 1:
            idx = np.argwhere(files==file)[0][0]
            scripts = self._flow[idx,1]
            if len(scripts) == 0: self._flow[idx,1] = script
            elif script in scripts: pass # do nothing if this input is already recorded
            else: self._flow[idx,1] += f', {script}'
        # if the file is not otherwise found in the workflow file
        elif len(np.argwhere(files==file)) == 0:
            self._flow = np.vstack((self._flow,
                np.array([file,'',script], dtype = '<U1024')))
        else: raise ValueError('Check workflow file, repetitive entries found')
        
        # save file
        np.savetxt(self._flow_file, self._flow, delimiter = '\t', fmt = '%s')