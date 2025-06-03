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
import argparse
import warnings

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
                 mode = DeprecationWarning('Use timeout instead of mode'), # deprecated, use timeout instead
                 lim = -1, # number of commands per file, default -1
                 arraysize = 100, # array size limit, default 2000 for CSD3 cluster, QOS max jobs 50
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
        
        import warnings
        if len(name) > 30:
            import secrets
            name = secrets.token_urlsafe(6).replace('-','_') # generate a random name
            warnings.warn(f'Job name too long, using random name {name}')
        self.name = name.replace('/','_') + '_0'
        self.debug = debug
        self.intr = intr
        self.parallel = parallel
        
        # SLURM config
        self.partition = partition
        self.timeout = timeout
        self.n_node = n_node
        self.n_task = n_task
        self.n_cpu = n_cpu
        arraysize = min(arraysize, int(500/n_cpu)) # QOS max CPU per user limit is 500
        self.arraysize = arraysize # the default array job size limit is 2000 for SLURM
        self.email = '--mail-type=ALL' if email else ''
        self.account = account
        
        # dependencies
        self.env = env
        self.wd = os.path.abspath(wd)
        if type(dependency) not in [list, tuple]: dependency = [dependency]
        self.dep = []
        for dep in dependency:
            if type(dep) in [int, array_submitter]: self.dep.append(dep)
        if type(modules) == type('a'): modules = [modules] # single string
        self.modules = modules
        
        # internal variables
        self._blank = True
        self._count = 1 # number of commands per file
        self._fileid = 0 # current file id
        self._nfiles = -1
        self._jobid = 0
        
        # status features
        self._submitted = False
        self._slurmid = []
        
        # limit of commands per file
        if lim > 0: self.lim = lim
        self._mode = 'long' if timeout > 15 else 'short' # long jobs are for > 15 min commands
        self.wallclock = 240 if timeout > 15 else 60 # timeout = 12 hrs max
        if wallclock > 0: self.wallclock = min(wallclock,720)
        self.lim = int(self.wallclock/timeout)
        self.lim = max(self.lim, 1) # at least one command per file
        self.lim *= self.parallel # all parallel processes have the same time limit
        
        # read command line args
        import __main__
        if 'args' in dir(__main__): self.config(**vars(__main__.args))
        print(f'Max {self.lim} commands per file, {self.arraysize} files per array job')

        # directories
        self.logdir = f'{log}/{self.name}/' # to prevent confusion with other array submissions
        if not os.path.isdir(self.logdir): os.mkdir(self.logdir)
        os.system(f'rm -rf {self.logdir}/*') # clear temp files from the previous run
        self.tmpdir = f'{tmpdir}/{self.name}/' # to prevent confusion with other array submissions
        if not os.path.isdir(self.tmpdir): os.mkdir(self.tmpdir)
        os.system(f'rm -rf {self.tmpdir}/*') # clear temp files from the previous run
    
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
            setattr(self, key, value)
        if 'jobname' in kwargs.keys() and kwargs['jobname'] != None and self._blank: 
            setattr(self, 'name', kwargs['jobname'])
        if 'wallclock' in kwargs.keys() and kwargs['wallclock'] != None:
            self.wallclock = min(self.wallclock,720)
            self.lim = int(self.wallclock/self.timeout)
            self.lim = max(self.lim, 1) # at least one command per file
            self.lim *= self.parallel # all parallel processes have the same time limit

    # a new file requires a shebang line, so this func resets the file
    def _newfile(self):
        fname = self.tmpdir+f'{self.name}_{self._fileid}.sh'
        _file = open(fname,'w')
        print('#!/bin/bash', file = _file) # shebang line
        if type(self.env) != type(None): # mamba environment
            print('source /home/yh464/.bashrc', file = _file)
            print(f'mamba activate {self.env}', file = _file)
        for mod in self.modules: # load modules
            print(f'module load {mod}', file = _file)
        print(f'cd {self.wd}', file = _file) # enforce working directory
        _file.close()
    
    def _newjob(self):
        # reset cmd count and file id
        self._count = 1
        self._fileid = 0
        self._jobid += 1
        self._nfiles = -1
        
        # reset job name and directories to avoid confusion
        self.name = self.name.replace(f'_{self._jobid-1}',f'_{self._jobid}')
        self.logdir = self.logdir.replace(f'_{self._jobid-1}',f'_{self._jobid}')
        if not os.path.isdir(self.logdir): os.mkdir(self.logdir)
        os.system(f'rm -rf {self.logdir}/*') # clear temp files from the previous run
        self.tmpdir = self.tmpdir.replace(f'_{self._jobid-1}',f'_{self._jobid}')
        if not os.path.isdir(self.tmpdir): os.mkdir(self.tmpdir)
        os.system(f'rm -rf {self.tmpdir}/*') # clear temp files from the previous run
        
    # adds a command
    def add(self,*cmd):
        '''
        Adds commands to the submitter utility
        Can pass any number of commands that need to be run in order
        '''
        # can append multiple commands at once if they need to be run successively
        # they are considered a single command and take care of the time limit!
        
        self._blank = False
        # preprocess command
        cmd = ' && '.join(cmd)
        if self.parallel > 1:
            cmd += ' &'
            if self._count % self.parallel == 0: cmd += ' wait'

        if self._mode == 'long':
            # first iteration: reset files
            if self._count == 1:
                self._nfiles += 1
                self._newfile()
            
            # append command
            fname = self.tmpdir+f'{self.name}_{self._fileid}.sh'
            _file = open(fname,'a')
            print(cmd, file = _file)
            _file.close() # I do not want 2000 file handles!
            
            # proceed to next file
            self._fileid += 1
            
            # if array size limit is reached, print a new command to the first file
            if self._fileid >= self.arraysize:
                if self._count >= self.lim:
                    self.submit() # because this job is 'full'
                    self._newjob()
                else:
                    self._count += 1
                    self._fileid = 0
        
        elif self._mode == 'short':
            # for 'short' jobs, submit in 15-min batches
            if self._count == 1:
                self._nfiles += 1
                self._newfile()
            
            fname = self.tmpdir+f'{self.name}_{self._fileid}.sh'
            _file = open(fname,'a')
            print(cmd, file = _file)
            _file.close() # I do not want 2000 file handles!
            
            # if the total time reaches the file size limit, start a new file
            if self._count >= self.lim:
                if self._fileid + 1 >= self.arraysize:
                    self.submit()
                    self._newjob()
                else:
                    self._count = 1
                    self._fileid += 1
            else:
                self._count += 1
    
    def _write_wrap(self):
        # master wrapper
        self._wrap_name = f'{self.tmpdir}/{self.name}_wrap.sh'       
        wrap = open(self._wrap_name,'w')
        
        # sbatch arguments
        print('#!/bin/bash',file = wrap)
        print(f'#SBATCH -N {self.n_node}', file = wrap)
        print(f'#SBATCH -n {self.n_task}', file = wrap)
        print(f'#SBATCH -c {self.n_cpu}', file = wrap)
        time = int(self.timeout * self._count / self.parallel)
        print(f'#SBATCH -t {time}', file = wrap)
        print(f'#SBATCH -p {self.partition}', file = wrap)
        print(f'#SBATCH -o {self.logdir}/{self.name}_%a.log', file = wrap)
        print(f'#SBATCH -e {self.logdir}/{self.name}_%a.err', file = wrap)
        print(f'#SBATCH --array=0-{self._nfiles}', file = wrap)
        dep_str = ''
        for dep in self.dep: 
            if type(dep) == int:
                dep_str += f':{dep}'
            elif type(dep) == array_submitter:
                if dep._blank: continue
                if not dep.submitted:
                    dep.submit()
                    print(f'Warning: {dep.name} is listed as a dependency and automatically submitted')
                for idx in dep._slurmid:
                    dep_str += f':{idx}'
        if len(dep_str) > 0: 
            print('#SBATCH --dependency afterok'+dep_str, file = wrap)
        
        if self.email: print('#SBATCH --mail-type=ALL', file = wrap)
        if self.account: print(f'#SBATCH -A {self.account}', file = wrap)
        print(f'bash {self.tmpdir}/{self.name}_'+'${SLURM_ARRAY_TASK_ID}.sh', file = wrap)
        wrap.close()
    
    def _print(self):
        '''
        Prints one sample command file and the sbatch command
        '''
        if self._nfiles < 0:
            print('No files to submit (!)')
            return
        
        import os
        os.system(f'cat {self.tmpdir}/{self.name}_0.sh')
        print(f'\n\nbash {self.tmpdir}/{self.name}_0.sh\n\n')
        self._write_wrap()
        
        # debug command also outputs the submit command
        time = int(self.timeout * self._count)
        email = '--mail-type=ALL' if self.email else ''
        account = f'-A {self.account}' if self.account else ''
        print(f'sbatch -N {self.n_node} -n {self.n_task} -c {self.n_cpu} '+
                  f'-t {time} -p {self.partition} {email} {account} '+
                  f'-o {self.logdir}/{self.name}_%a.log -e {self.logdir}/{self.name}_%a.err'+ # %a = array index
                  f' --array=0-{self._nfiles} {self._wrap_name}') 
    
    def _splash(self, jobid):
        msg = []
        msg.append('#' * 100)
        msg.append('Following job has been submitted to SLURM:')
        msg.append(f'    Name:      {self.name}')
        msg.append(f'    Path:      {self.tmpdir}')
        msg.append(f'    Partition: {self.partition}')
        msg.append(f'    Timeout:   {int(self.timeout*self._count/self.parallel)} minutes')
        msg.append(f'    CPUs:      {self.n_cpu}')
        msg.append(f'    # files:   {self._nfiles + 1}')
        msg.append(f'    Parallel:  {self.parallel}')
        msg.append(f'    Job ID:    {jobid}')
        msg.append('#' * 100)
        if not self._blank: print('\n'.join(msg))

    def submit(self):
        '''
        Submits all commands to the cluster (SLURM manager)
        '''
        # if debug mode is on, debug instead
        if self._blank: return
        if self.debug: self._print(); return
        if self.intr: import os; os.system(f'for x in {self.tmpdir}/*.sh; do bash $x; done'); return
        
        self._write_wrap()
        time = int(self.timeout) * self._count
        email = '--mail-type=ALL' if self.email else ''
        account = f'-A {self.account}' if self.account else ''
        
        from subprocess import check_output
        msg = check_output(f'sbatch -N {self.n_node} -n {self.n_task} -c {self.n_cpu} '+
                  f'-t {time} -p {self.partition} {email} {account} '+
                  f'-o {self.logdir}/{self.name}_%a.log -e {self.logdir}/{self.name}_%a.err'+ # %a = array index
                  f' --array=0-{self._nfiles} {self._wrap_name}', shell = True
                  ).decode().strip()
        jobid = int(msg.split()[-1]) # raises an error if sbatch fails
        self._splash(jobid)
        self._slurmid.append(jobid)
        self.submitted = True
        return jobid

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
        slurm.add_argument('--timeout', help = 'timeout in minutes')
        slurm.add_argument('--wallclock', help = 'total time limit per file in minutes')
        slurm.add_argument('--parallel', help = 'number of parallel processes')
        slurm.add_argument('--debug', help = 'debug mode', default = False, action = 'store_true')
        slurm.add_argument('--dep', help = 'dependencies', default = [], nargs = '*')
        slurm.add_argument('--intr', help = 'interactive mode', default = False, action = 'store_true')
        