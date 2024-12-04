'''
This script is a utility to log the command history
'''

import os
import datetime as dt
import sys

default_dir = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/cmdhistory/'
default_filename = 'current_cmdhistory.txt'

def newfile():
    now = dt.datetime.now(dt.timezone.utc)
    date = now.date().isoformat()
    if os.path.isfile(default_dir+default_filename):
        os.system(f'mv {default_dir}{default_filename} {default_dir}{date}_cmdhistory.txt')
    f = open(default_dir+default_filename, 'w')
    now = now.strftime('%Y-%m-%d %H:%M:%S')
    print(f'Command history from UTC {now}: \n', file = f)
    f.close()

def log():
    # prints the current command to the cmdhistory file, verbatim
    now = dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%d %H:%M:%S')
    if not os.path.isfile(default_dir+default_filename):
        newfile()
    f = open(default_dir+default_filename, 'a')
    print(now, file = f)
    print('python '+' '.join(sys.argv),file = f)
    print('WD = '+os.getcwd(), file = f)
    print(file = f)
    f.close()
    # automatically sets up a new file when the current file reaches a threshold size
    # 100000 bytes (can be changed)
    if os.path.getsize(default_dir+default_filename) > 100000:
        newfile()
    
def errlog():
    # prints the current error to the cmdhistory file
    f = open(default_dir+default_filename, 'a')
    print('Above command raised an error', file = f)
    for x in sys.exc_info(): # exception info
        print(repr(x),file = f)
    print(file = f)
    raise
    # requires the following design:
    # try: main(args)
    # except: cmdhistory.errlog()