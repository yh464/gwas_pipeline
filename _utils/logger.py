#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2024-11-25

This is a general utility to print splash screens with input command-line
arguments and output a log file
'''
import os, sys, datetime, inspect
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
    def __init__(self, fname = None, echo = True, **kwargs):
        self.file = open(fname, 'w') if fname is not None else sys.stdout
        self.echo = echo
        
    def log(self, msg, warning = False):
        now = datetime.datetime.now().isoformat(sep = ' ')
        calling_file = os.path.basename(inspect.stack()[1].filename).replace('.py','')
        warning_str = '| WARNING ' if warning else ''
        msg = f'[ {now} | {calling_file} {warning_str}] {msg}'
        print(msg, file = self.file)
        if self.echo and self.file != sys.stdout: print(msg)
    
    def splash(self, args):
        msg = splash(args, silent = True)
        self.log(msg)