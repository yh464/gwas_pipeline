#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2024-11-25
Version 2: 2025-12-26

This is a general utility to print splash screens with input command-line
arguments and output a log file
'''
import os, sys, datetime, inspect, tracemalloc, time
def splash(args, silent = False):
    msg = []
    msg.append('=' * 100)
    msg.append('Calling script:')
    msg.append(f'    {sys.argv[0]}')
    slurm_args = ['jobname','name','debug','partition', 'timeout', 'n_node', 'n_task','n_cpu', 'n_gpu',
            'arraysize','email','account','env','wd','dep','modules','logdir',
            'tmpdir','lim','intr','wallclock','parallel']
    msg.append('Input options:')
    v = vars(args)
    for var in v:
        if var in slurm_args: continue
        val = v[var]
        msg.append(f'    {var!s:25}{val!s}')
    
    slurm = []
    slurm.append('Slurm management options:')
    for var in v:
        if var not in slurm_args: continue
        val = v[var]
        if val == None: val = '(default)'
        slurm.append(f'    {var!s:25}{val!s}')
    if len(slurm) > 1: msg += slurm
    msg.append('=' * 100)
    if not silent: print('\n'.join(msg))
    return '\n'.join(msg)

class logger():
    def __init__(self, fname = None, echo = True, **kwargs):
        self.file = open(fname, 'w') if fname is not None else sys.stdout
        self.echo = echo
        self.start_time = time.perf_counter()
        self.cpu_time = time.process_time()
        tracemalloc.start()
        
    def log(self, msg, warning = False, error = False, info = False, calling_file = None):
        now = datetime.datetime.now().isoformat(sep = ' ')
        if calling_file is None:
            calling_file = os.path.basename(inspect.stack()[2].filename).replace('.py','')
        if error: warning_str = '| ERROR '
        elif warning: warning_str = '| WARNING '
        elif info: warning_str = '| INFO '
        else: warning_str = ''
        msg = f'[ {now} | {calling_file} {warning_str}] {msg}'
        print(msg, file = self.file)
        if self.echo and self.file != sys.stdout: print(msg)
    
    def warn(self, msg, calling_file = None): self.log(msg, warning = True, calling_file = calling_file)
    def error(self, msg, calling_file = None): self.log(msg, error = True, calling_file = calling_file)
    def info(self, msg, calling_file = None): self.log(msg, info = True, calling_file = calling_file)

    def system(self, cmd, calling_file = None):
        out = os.system(cmd)
        self.log(cmd, calling_file = calling_file)
        self.log(f'Exit code: {out}', calling_file = calling_file)
        return out

    def check_memory(self, calling_file = None):
        mem = tracemalloc.get_traced_memory()[0] / 1024 ** 2 
        self.log(f'Current memory usage: {mem:.2f} MB', calling_file = calling_file) 
        return mem

    def splash(self, args):
        msg = splash(args, silent = True)
        print(msg, file = self.file)
        if self.echo and self.file != sys.stdout: print(msg)

    def __del__(self):
        self.log('Analysis finished')
        self.log(f'    Peak memory usage: {tracemalloc.get_traced_memory()[1] / 1024 ** 2:.2f} MB')
        self.log(f'    Total time: {time.perf_counter() - self.start_time:.2f} seconds')
        self.log(f'    CPU time: {time.process_time() - self.cpu_time:.2f} seconds')
        self.log(f'    CPU usage: {(time.process_time() - self.cpu_time) / (time.perf_counter() - self.start_time) * 100:.2f}%')
        tracemalloc.stop()