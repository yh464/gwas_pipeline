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
    msg.append('Input options:')
    v = vars(args)
    for var in v:
        val = v[var]
        msg.append(f'    {var}' + ' ' * (10-len(var)) + str(val))
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