# -*- coding: utf-8 -*-
#!/usr/env/bin python3

'''
This utility manages the nomenclature of file names / paths
'''

import os
import numpy as np
import re

class project():
    def __init__(self,
                 _dir = '../path/', # uses relative path, so defaults to the relative path to the wd
                 _dict = 'dict.txt', # nomenclature
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