'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2024-08-22
Version 2: 2025-09-25

This script is a general scaffold to submit jobs in batches
NB each command is executed separately in the bash script, so
if slurm returns 'FAILED', some steps may still run normally
CHECK LOG
'''

from logging import warning
import os
import argparse
import re
import warnings
import math
from hashlib import sha256

class array_submitter():
    '''
    Attributes of an array submitter
    Required:
        Name (must be unique, will overwrite other submitters)
        Timeout (per command or group of commands)
    Optional:
        n_node, n_task, n_cpu (resource allocation)
        log (log file directory)
        tmpdir (stores batch scripts)
        lim (limit of commands per file)
        arraysize (limit of number of files before starting a new array job)
        email: True/False
        mode: 'short' or 'long' jobs
        env (mamba environment)
        modules (module load *)
        dependency (can be another submitter or an array job ID)
        account
        wd (working directory)
        debug: True/False
    '''
    # NB update default resource settings by checking sacctmgr list QOS
    def __init__(self,
                 name, # name of project
                 timeout, # time limit per command, minutes; will raise an error if not given
                 partition = 'icelake-himem',
                 n_node = 1,
                 n_task = 1,
                 n_cpu = 1,
                 log = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/logs',
                 tmpdir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/temp',
                 parallel = 1, # number of parallel processes, useful for small jobs that need <1 CPU
                 lim = warnings.warn(DeprecationWarning('Command limit will be automatically determined')), # deprecated
                 arraysize = 200, # array size limit, default 2000 for CSD3 cluster, QOS max jobs 500
                 email = True,
                 wallclock = -1, # total time limit per file, default 240 minutes
                 env = 'wd', # default working environment
                 modules = [], # modules to load
                 dependency = [], # dependent jobs
                 account = None,
                 wd = '.',
                 debug = False,
                 intr = False # tries to run interactively in series (CAUTION WITH THIS OPTION)
                 ):

        # current SLURM limit: 50 jobs, 450 CPUs; adjust for busy cluster
        cpu_avail = {
            'desktop': 4,
            'sapphire': 1,
            'icelake-himem': 1,
            'icelake': 1,
            'cclake': 1,
            'cclake-himem': 1
        }

        self._full_name = name
        if len(name) > 30:
            name = sha256(name.encode()).hexdigest()[:6] # truncate to 6S characters
            warnings.warn(f'Job name too long, using random name {name}')
        self.name = '_' + name.replace('/','_') + '_0'
        self.debug = debug
        self.intr = intr
        self.parallel = parallel
        
        # SLURM config
        self.partition = partition
        self.timeout = timeout
        self.n_node = n_node
        self.n_task = n_task
        self.n_cpu = n_cpu
        if n_cpu < 20: arraysize = min(arraysize, int(500/n_cpu)) # QOS max CPU per user limit is 500
        self.arraysize = arraysize # the default array job size limit is 2000 for SLURM
        self.email = '--mail-type=ALL' if email else ''
        self.account = account
        
        # dependencies
        self.env = env
        self.wd = os.path.abspath(wd)
        if type(dependency) not in [list, tuple]: dependency = [dependency]
        self.dep = []
        for dep in dependency:
            if type(dep) == array_submitter: self.dep.append(dep)
            else: self.dep.append(int(dep))
        if type(modules) == type('a'): modules = [modules] # single string
        self.modules = modules
        
        # limit of commands per file and wallclock limit
        self.wallclock = max(timeout, 240) if timeout > 15 else 60
        if wallclock > 0: self.wallclock = wallclock
        if account != None and account.find('sl2') >= 0: self.wallclock = min(self.wallclock, 2160)
        else: self.wallclock = min(self.wallclock, 720) 

        # read command line args before specifying limit of commands per file
        import __main__
        if 'args' in dir(__main__): self.config(**vars(__main__.args))

        # adjust parallel processes based on available CPUs
        if self.n_cpu < cpu_avail[self.partition]:
            self.parallel = math.floor(self.parallel * cpu_avail[self.partition]/self.n_cpu)
            self.n_cpu = cpu_avail[self.partition]

        # number of *parallel batches* of commands per file
        self.lim = int(self.wallclock/timeout)
        self.lim = max(self.lim, 1) # at least one command per file
        print(f'Max {self.lim} batches * {self.parallel} commands per file, {self.arraysize} files per array job')

        # directories
        self.logdir = f'{log}/{self.name}' # to prevent confusion with other array submissions
        if not os.path.isdir(self.logdir): os.mkdir(self.logdir)
        os.system(f'rm -rf {self.logdir}/*') # clear temp files from the previous run
        self.tmpdir = f'{tmpdir}/{self.name}' # to prevent confusion with other array submissions
        if not os.path.isdir(self.tmpdir): os.mkdir(self.tmpdir)
        os.system(f'rm -rf {self.tmpdir}/*') # clear temp files from the previous run

        # commands are staged up to an array size limit before a new job array is initialised
        self.array_cmd_limit = self.arraysize * self.lim * self.parallel
        self._staged_cmd = []
        self._blank = True
        self._count = 0
        self._nfiles = 0
        self._jobid = 0
        
        # SLURM status properties
        self.submitted = False
        self._slurmid = []
    
    # change settings
    def config(self, **kwargs):
        valid_keys = [
            'name', 'debug','partition', 'timeout','wallclock','n_node', 'n_task','n_cpu',
            'arraysize','email','account','env','wd','dep','modules','logdir',
            'tmpdir','lim','intr','dependency', 'parallel'
            ]
        numeric_keys = ['n_node', 'n_task', 'n_cpu', 'timeout', 'arraysize', 'lim', 'wallclock','parallel']
        for key, value in kwargs.items():
            if key in numeric_keys and value != None and value.startswith('x'):
                # if the value starts with 'x', it is a multiplier
                value = float(value[1:]) * getattr(self, key)
            if not key in valid_keys: continue
            if value == None: continue # skip None values
            if key in numeric_keys: value = int(value) # convert to int
            if key in ['dep','dependency']: value = self.dep + value # dependencies are appended and not overridden
            setattr(self, key, value)
        if 'jobname' in kwargs.keys() and kwargs['jobname'] != None and len(self._staged_cmd) == 0: 
            setattr(self, 'name', kwargs['jobname'])
        if 'wallclock' in kwargs.keys() and kwargs['wallclock'] != None:
            self.wallclock = min(self.wallclock,720)
        self.lim = int(self.wallclock/self.timeout)
        self.lim = max(self.lim, 1) # at least one command per file
        self.lim *= self.parallel # all parallel processes have the same time limit

    # adds a command
    def add(self,*cmd):
        '''
        Adds commands to the submitter utility
        Can pass any number of commands that need to be run in the specific order
        '''
        # can append multiple commands at once if they need to be run successively
        # they are considered a single command and take care of the time limit!
        if self.submitted: raise ValueError('Cannot add commands after submission')
        if len(cmd) == 0: return
        elif len(cmd) > 1: cmd = '(' + ' && '.join(cmd) + ')'; self._staged_cmd.append(cmd)
        else: cmd = cmd[0]; self._staged_cmd.append(cmd)
        self._blank = False

    # initialises a new job array
    def _newjob(self):
        # reset cmd count and file id
        self._count = 0
        self._nfiles = 0
        self._jobid += 1
        self._wrap_name = ''
        
        # reset job name and directories to avoid confusion
        self.name = re.sub(f'_{self._jobid-1}$',f'_{self._jobid}', self.name)
        self.logdir = os.path.realpath(f'{self.logdir}/..')+f'/{self.name}'
        if not os.path.isdir(self.logdir): os.mkdir(self.logdir)
        os.system(f'rm -rf {self.logdir}/*') # clear temp files from the previous run
        self.tmpdir = os.path.realpath(f'{self.tmpdir}/..')+f'/{self.name}'
        if not os.path.isdir(self.tmpdir): os.mkdir(self.tmpdir)
        os.system(f'rm -rf {self.tmpdir}/*') # clear temp files from the previous run

    # dumps staged commands to a new file
    def _write_file(self, fileid, cmds = []):
        fname = f'{self.tmpdir}/{self.name}_{fileid}.sh'
        _file = open(fname,'w')
        print('#!/bin/bash', file = _file) # shebang line
        if type(self.env) != type(None): # mamba environment
            print('source ~/.bashrc', file = _file)
            print(f'mamba activate {self.env}', file = _file)
        for mod in self.modules: # load modules
            print(f'module load {mod}', file = _file)
        print(f'cd {self.wd}', file = _file) # enforce working directory
        for cmd in cmds: print(cmd, file = _file) # print commands
        _file.close()
    
    # generates the dependency string
    def _write_dep_str(self):
        dep_str = []
        for dep in self.dep: 
            if type(dep) in [int, str]:
                dep_str.append(f'afterok:{dep}')
            else:
                if dep._blank: 
                    print(f'Warning: {dep.name} has no commands to run, skipping dependency')
                    continue
                if not dep.submitted:
                    dep.submit()
                    print(f'Warning: {dep.name} is listed as a dependency and automatically submitted')
                for idx in dep._slurmid:
                    dep_str.append(f'afterok:{idx}')
                if len(dep._slurmid) == 0:
                    dep_str.append('<please run again removing the --debug flag>')
        if len(dep_str) > 0: 
            return '-d '+','.join(dep_str)
        else: return ''

    # writes the wrapper to be submitted using sbatch
    def _write_wrap(self):
        # master wrapper
        self._wrap_name = f'{self.tmpdir}/{self.name}_wrap.sh'       
        wrap = open(self._wrap_name,'w')
        if len(self._staged_cmd) < self.parallel: n_cpu = math.ceil(self.n_cpu / self.parallel * len(self._staged_cmd))
        else: n_cpu = self.n_cpu

        # sbatch arguments
        print('#!/bin/bash',file = wrap)
        print(f'#SBATCH -N {self.n_node}', file = wrap)
        print(f'#SBATCH -n {self.n_task}', file = wrap)
        print(f'#SBATCH -c {n_cpu}', file = wrap)
        time = self.timeout * self._count
        print(f'#SBATCH -t {int(time)}', file = wrap)
        print(f'#SBATCH -p {self.partition}', file = wrap)
        print(f'#SBATCH -o {self.logdir}/{self.name}_%a.log', file = wrap)
        print(f'#SBATCH -e {self.logdir}/{self.name}_%a.err', file = wrap)
        print(f'#SBATCH --array=0-{self._nfiles-1}', file = wrap)
        print('#SBATCH '+ self._write_dep_str(), file = wrap)
        
        if self.email: print('#SBATCH --mail-type=ALL', file = wrap)
        if self.account: print(f'#SBATCH -A {self.account}', file = wrap)
        print(f'bash {self.tmpdir}/{self.name}_'+'${SLURM_ARRAY_TASK_ID}.sh', file = wrap)
        wrap.close()

    # dumps staged commands to a job array
    def _dump(self):
        if len(self._staged_cmd) == 0: return
        cmds_to_dump = self._staged_cmd[:min(len(self._staged_cmd), self.array_cmd_limit)]

        # organise cmds into parallel batches
        if self.parallel > 1:
            batches = []
            for i in range(0, len(cmds_to_dump), self.parallel):
                batch = cmds_to_dump[i:min(i+self.parallel, len(cmds_to_dump))]
                batch = ' & '.join(batch) + '; for job in $(jobs -p); do wait $job; done' # run in parallel
                batches.append(batch)
            cmds_to_dump = batches; del batches

        self._nfiles = min(len(cmds_to_dump), self.arraysize)
        cmds_to_dump = [cmds_to_dump[i::self._nfiles] for i in range(self._nfiles)]
        for fileid, cmds in enumerate(cmds_to_dump):
            self._write_file(fileid, cmds)
        self._count = len(cmds_to_dump[0]) # number of parallel batches per file
        self._write_wrap()

        if len(self._staged_cmd) > self.array_cmd_limit:
            self._submit_single()
            self._staged_cmd = self._staged_cmd[self.array_cmd_limit:]
            self._newjob(); self._dump()
    
    # for debugging
    def _print(self):
        '''
        Prints one sample command file and the sbatch command
        '''
        if self._nfiles == 0:
            print('No files to submit (!)')
            return
        
        print('\n' + '#' * 100 + '\n')

        print(f'Job name: {self._full_name}')

        import os
        os.system(f'cat {self.tmpdir}/{self.name}_0.sh')
        print(f'\n\nbash {self.tmpdir}/{self.name}_0.sh\n\n')
        
        dep_str = self._write_dep_str() # this will submit dependencies which will be printed first in debug mode

        # debug command also outputs the submit command
        time = self.timeout * self._count
        email = '--mail-type=ALL' if self.email else ''
        account = f'-A {self.account}' if self.account else ''
        print(f'sbatch -N {self.n_node} -n {self.n_task} -c {self.n_cpu} '+
                  f'-t {int(time)} -p {self.partition} {email} {account} '+
                  f'-o {self.logdir}/{self.name}_%a.log -e {self.logdir}/{self.name}_%a.err'+ # %a = array index
                  f' --array=0-{self._nfiles-1} {dep_str} {self._wrap_name}') 
        print('\n' + '#' * 100 + '\n')
    
    # splash screen
    def _splash(self, jobid):
        if len(self._staged_cmd) < self.parallel: n_cpu = math.ceil(self.n_cpu / self.parallel * len(self._staged_cmd))
        else: n_cpu = self.n_cpu
        msg = []
        msg.append('#' * 100)
        msg.append('Following job has been submitted to SLURM:')
        msg.append(f'    Name:       {self._full_name}')
        msg.append(f'    Path:       {self.tmpdir}')
        msg.append(f'    Partition:  {self.partition}')
        msg.append(f'    Timeout:    {self.timeout * math.ceil(self._count / self.parallel)} minutes')
        msg.append(f'    CPUs:       {n_cpu}')
        msg.append(f'    # files:    {self._nfiles}')
        msg.append(f'    Parallel:   {self.parallel}')
        msg.append(f'    Dependency: {re.sub("^-d ","", self._write_dep_str()) if self._write_dep_str() != "" else "None"}')
        msg.append(f'    Job ID:     {jobid}')
        msg.append('#' * 100)
        if not self._blank: print('\n'.join(msg))

    # submits a single job array
    def _submit_single(self):
        if self.debug: self._print(); self.submitted = True; return # debug mode -> print only
        if self.intr: os.system(f'for x in {self.tmpdir}/*.sh; do bash $x; done'); return
        time = self.timeout * math.ceil(self._count / self.parallel)
        time = min(time, 720)
        email = '--mail-type=ALL' if self.email else ''
        account = f'-A {self.account}' if self.account else ''
        
        from subprocess import check_output
        msg = check_output(f'sbatch -N {self.n_node} -n {self.n_task} -c {self.n_cpu} '+
                  f'-t {int(time)} -p {self.partition} {email} {account} '+
                  f'-o {self.logdir}/{self.name}_%a.log -e {self.logdir}/{self.name}_%a.err'+ # %a = array index
                  f' --array=0-{self._nfiles-1} {self._write_dep_str()} {self._wrap_name}', shell = True
                  ).decode().strip()
        jobid = int(msg.split()[-1]) # raises an error if sbatch fails
        self._splash(jobid)
        self._slurmid.append(jobid)
        self.submitted = True

    # master wrapper for job submission
    def submit(self):
        '''
        Submits all commands to the cluster (SLURM manager)
        '''
        if self._blank: return
        if self.submitted: warnings.warn(f'Job {self.name} already submitted'); return
        self._dump()
        self._submit_single()
        self._staged_cmd = []

class slurm_parser(argparse.ArgumentParser):
    '''
    An argparse.ArgumentParser with default SLURM config options
    '''
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.parser_config()

    def parser_config(self):
        slurm = self.add_argument_group('SLURM configuration, enter numbers to override default resource allocation,\n'+
            'enter x2, etc. to multiply the default values, leave blank to use defaults')
        slurm.add_argument('--jobname', help = 'Manually specify job name')
        slurm.add_argument('--partition', help = 'partition')
        slurm.add_argument('--account', help = 'account to charge')
        slurm.add_argument('--n_cpu', help = 'number of CPUs per task')
        slurm.add_argument('--n_node', help = 'number of nodes needed')
        slurm.add_argument('--n_task', help = 'number of tasks per job')
        slurm.add_argument('--arraysize', help = 'number of files per array job')
        slurm.add_argument('--timeout', help = 'timeout in minutes')
        slurm.add_argument('--wallclock', help = 'total time limit per file in minutes')
        slurm.add_argument('--parallel', help = 'number of parallel processes')
        slurm.add_argument('--debug', help = 'debug mode', default = False, action = 'store_true')
        slurm.add_argument('--dep', help = 'dependencies', default = [], nargs = '*')
        slurm.add_argument('--intr', help = 'interactive mode', default = False, action = 'store_true')
        