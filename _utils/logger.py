#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2024-11-25

This is a general utility to print splash screens with input command-line
arguments and output a log file
'''

def splash(args, silent = False):
    import sys
    msg = []
    msg.append('=' * 100)
    msg.append('Calling script:')
    msg.append(f'    {sys.argv[0]}')
    slurm_args = ['name','debug','partition', 'timeout', 'n_node', 'n_task','n_cpu',
            'arraysize','email','account','env','wd','dep','modules','logdir',
            'tmpdir','lim','intr']
    msg.append('Input options:')
    v = vars(args)
    for var in v:
        if var in slurm_args: continue
        val = v[var]
        msg.append(f'    {var}' + ' ' * (15-len(var)) + str(val))
    
    slurm = []
    slurm.append('Slurm management options:')
    for var in v:
        if var not in slurm_args: continue
        val = v[var]
        slurm.append(f'    {var}' + ' ' * (15-len(var)) + str(val))
    if len(slurm) > 1: msg += slurm
    msg.append('=' * 100)
    if not silent: print('\n'.join(msg))
    return '\n'.join(msg)

class logger():
    def __init__(self, fname, silent = False):
        self.log = open(fname, 'w')
        self.silent = silent
        
    def log(self, msg):
        print(msg, file = self.log)
        if not self.silent: print(msg)
    
    def splash(self, args):
        msg = splash(args, silent = True)
        self.log(msg)