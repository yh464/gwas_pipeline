#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2024-11-25

This is a general utility to print splash screens with input command-line
arguments and output a log file
'''
import os, sys
def splash(args, silent = False):
    msg = []
    msg.append('=' * 100)
    msg.append('Calling script:')
    msg.append(f'    {sys.argv[0]}')
    slurm_args = ['jobname','name','debug','partition', 'timeout', 'n_node', 'n_task','n_cpu',
            'arraysize','email','account','env','wd','dep','modules','logdir',
            'tmpdir','lim','intr','wallclock','parallel']
    msg.append('Input options:')
    v = vars(args)
    for var in v:
        if var in slurm_args: continue
        val = v[var]
        msg.append(f'    {var!s:15}{val!s}')
    
    slurm = []
    slurm.append('Slurm management options:')
    for var in v:
        if var not in slurm_args: continue
        val = v[var]
        if val == None: val = '(default)'
        slurm.append(f'    {var!s:15}{val!s}')
    if len(slurm) > 1: msg += slurm
    msg.append('=' * 100)
    if not silent: print('\n'.join(msg))
    return '\n'.join(msg)

class logger():
    def __init__(self, fname = None, **kwargs):
        self.file = open(fname, 'w') if fname is not None else sys.stdout
        if 'silent' in kwargs.keys(): self.silent = kwargs['silent']
        else: self.silent = False
        
    def log(self, msg):
        print(msg, file = self.file)
        if not self.silent: print(msg)
    
    def splash(self, args):
        msg = splash(args, silent = True)
        self.log(msg)