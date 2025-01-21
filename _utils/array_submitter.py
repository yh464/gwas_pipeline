'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2024-08-22

This script is a general scaffold to submit jobs in batches
NB each command is executed separately in the bash script, so
if slurm returns 'FAILED', some steps may still run normally
CHECK LOG
'''

import os

class array_submitter():
    def __init__(self,
                 name, # name of project
                 timeout, # time limit per command, minutes; will raise an error if not given
                 partition = 'icelake-himem',
                 n_node = 1,
                 n_task = 1,
                 n_cpu = 1,
                 log = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/logs',
                 tmpdir = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/temp',
                 lim = -1, # number of commands per file, default -1
                 arraysize = 500, # array size limit, default 2000 for CSD3 cluster, QOS max CPU per user limit 500
                 email = True,
                 mode = None,
                 env = 'wd', # default working environment
                 modules = [], # modules to load
                 account = None,
                 wd = '.',
                 debug = False,
                 ):
        
        self.name = name + '_0'
        self.partition = partition
        self.timeout = timeout
        self.n_node = n_node
        self.n_task = n_task
        self.n_cpu = n_cpu
        self.env = env
        self.wd = os.path.abspath(wd)
        arraysize = min(arraysize, int(500/n_cpu)) # QOS max CPU per user limit is 500
        
        self.logdir = f'{log}/{self.name}/' # to prevent confusion with other array submissions
        if not os.path.isdir(self.logdir): os.mkdir(self.logdir)
        os.system(f'rm -rf {self.logdir}/*') # clear temp files from the previous run
        self.tmpdir = f'{tmpdir}/{self.name}/' # to prevent confusion with other array submissions
        if not os.path.isdir(self.tmpdir): os.mkdir(self.tmpdir)
        os.system(f'rm -rf {self.tmpdir}/*') # clear temp files from the previous run
        
        if type(modules) == type('a'): modules = [modules] # single string
        self.modules = modules
        
        # internal variables
        self._email = '--mail-type=ALL' if email else ''
        self._account = f'-A {account}' if type(account) != type(None) else ''
        self._arraysize = arraysize # the default array job size limit is 2000 for SLURM
        self._count = 1 # number of commands per file
        self._fileid = 0 # current file id
        self._jobid = 0
        self._mode = 'long' if timeout > 15 else 'short'
        if type(mode) != type(None): self._mode = mode # override mode option if given
        self._debug = debug
        
        # limit per file
        max_time = 720 if self._mode == 'long' else 15 # timeout = 12 hrs max
        if lim != -1: max_time = 720 # so that time limit can be manually specified
        import math
        if lim == -1:
            self.lim = math.floor(max_time/timeout)
        else:
            self.lim = min((lim, math.floor(max_time/timeout)))
        del math, max_time
        
    # a new file requires a shebang line, so this func resets the file
    def _newfile(self):
        fname = self.tmpdir+f'{self.name}_{self._fileid}.sh'
        _file = open(fname,'w')
        print('#!/bin/bash', file = _file) # shebang line
        if type(self.env) != type(None): # conda environment
            print('source /home/yh464/.bashrc', file = _file)
            print(f'conda activate {self.env}', file = _file)
        for mod in self.modules: # load modules
            print(f'module load {mod}', file = _file)
        print(f'cd {self.wd}', file = _file) # enforce working directory
        _file.close()
    
    def _newjob(self):
        # reset cmd count and file id
        self._count = 1
        self._fileid = 0
        self._jobid += 1
        
        # reset job name and directories to avoid confusion
        self.name = self.name.replace(f'_{self._jobid-1}',f'_{self._jobid}')
        self.logdir = self.logdir.replace(f'_{self._jobid-1}',f'_{self._jobid}')
        if not os.path.isdir(self.logdir): os.mkdir(self.logdir)
        os.system(f'rm -rf {self.logdir}/*') # clear temp files from the previous run
        self.tmpdir = self.tmpdir.replace(f'_{self._jobid-1}',f'_{self._jobid}')
        if not os.path.isdir(self.tmpdir): os.mkdir(self.tmpdir)
        os.system(f'rm -rf {self.tmpdir}/*') # clear temp files from the previous run
        
    # adds a command
    def add(self,cmd):
        if self._mode == 'long':
            # first iteration: reset files
            if self._count == 1:
                self._newfile()
            
            # append command
            fname = self.tmpdir+f'{self.name}_{self._fileid}.sh'
            _file = open(fname,'a')
            print(cmd, file = _file)
            _file.close() # I do not want 2000 file handles!
            
            # proceed to next file
            self._fileid += 1
            
            # if array size limit is reached, print a new command to the first file
            if self._fileid >= self._arraysize:
                if self._count >= self.lim:
                    self.submit() # because this job is 'full'
                    self._newjob()
                else:
                    self._count += 1
                    self._fileid = 0
        
        elif self._mode == 'short':
            # for 'short' jobs, submit in 15-min batches
            if self._count == 1:
                self._newfile()
            
            fname = self.tmpdir+f'{self.name}_{self._fileid}.sh'
            _file = open(fname,'a')
            print(cmd, file = _file)
            _file.close() # I do not want 2000 file handles!
            
            # if the total time reaches the file size limit, start a new file
            if self._count >= self.lim:
                if self._fileid >= self._arraysize:
                    self.submit()
                    self._newjob()
                else:
                    self._count = 1
                    self._fileid += 1
            else:
                self._count += 1
    
    def debug(self):
        import os
        from fnmatch import fnmatch
        nfiles = 0
        for x in os.listdir(self.tmpdir):
            if fnmatch(x, f'{self.name}_*.sh'):
                n = x.replace(self.name,'').replace('.sh','').replace('_','')
                n = int(n)
                if n > nfiles: nfiles = n
        for i in range(min((int(5/self._count),nfiles))):
            os.system(f'cat {self.tmpdir}/{self.name}_{i}.sh')
        
        # master wrapper
        wrap_name = f'{self.tmpdir}/{self.name}_wrap.sh'
        wrap = open(wrap_name,'w')
        print('#!/bin/bash',file = wrap)
        print(f'bash {self.tmpdir}/{self.name}_'+'${SLURM_ARRAY_TASK_ID}.sh', file = wrap)
        wrap.close()
        
        # debug command also outputs the submit command
        time = self.timeout * self._count
        
        if nfiles < 0: return
        print(f'sbatch -N {self.n_node} -n {self.n_task} -c {self.n_cpu} '+
                  f'-t {time} -p {self.partition} {self._email} {self._account} '+
                  f'-o {self.logdir}/{self.name}_%a.log -e {self.logdir}/{self.name}_%a.err'+ # %a = array index
                  f' --array=0-{nfiles} {wrap_name}') 
    
    def submit(self):
        # if debug mode is on, debug instead
        if self._debug:
            self.debug()
            return
        
        # master wrapper
        wrap_name = f'{self.tmpdir}/{self.name}_wrap.sh'
        wrap = open(wrap_name,'w')
        print('#!/bin/bash',file = wrap)
        print(f'bash {self.tmpdir}/{self.name}_'+'${SLURM_ARRAY_TASK_ID}.sh', file = wrap)
        wrap.close()
        
        # scans directory for number of files (last sanity check)
        import os
        from fnmatch import fnmatch
        nfiles = 0
        for x in os.listdir(self.tmpdir):
            if fnmatch(x, f'{self.name}_*.sh'):
                n = x.replace(self.name,'').replace('.sh','').replace('_','')
                n = int(n)
                if n > nfiles: nfiles = n
        time = self.timeout * self._count
        
        if nfiles < 0: return
        os.system(f'sbatch -N {self.n_node} -n {self.n_task} -c {self.n_cpu} '+
                  f'-t {time} -p {self.partition} {self._email} {self._account} '+
                  f'-o {self.logdir}/{self.name}_%a.log -e {self.logdir}/{self.name}_%a.err'+ # %a = array index
                  f' --array=0-{nfiles} {wrap_name}') 