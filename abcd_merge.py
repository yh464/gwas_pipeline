# This script moves all files for one subject in ABCD
# from the datalad project to the existing Data_Imaging directory

import argparse
parser = argparse.ArgumentParser(description='Moves all files for single subject')
parser.add_argument('subj', help = 'Subject ID')
parser.add_argument('-i', dest = '_in', help = 'input directory',
    default = '/rds/project/rds-CeXlNYOYMxw/userdata/yh464/ABCD_fMRIprep/fmriprep')
parser.add_argument('-o', dest = 'out', help = 'target directory',
    default = '/rds/project/rds-CeXlNYOYMxw/Data_Imaging')
args=parser.parse_args()

import os
from fnmatch import fnmatch
import re

tmpdir = '/rds/project/rds-CeXlNYOYMxw/userdata/yh464/temp'
if not os.path.isdir(tmpdir): os.mkdir(tmpdir)

os.chdir(args._in)

# datalad project only contains symbolic links, so we need to extract contents first
os.system(f'datalad export-archive -d {args.subj} -c "" --missing-content ignore {tmpdir}/{args.subj}')
os.chdir(tmpdir)
os.system(f'tar -xvf {args.subj}.tar')
os.chdir(args.subj)

for d in os.listdir():
    if not fnmatch(d, 'ses-*'): # we only care about imaging data
        continue
    target_d = re.sub(r'Arm.*','',d)
    target_d = target_d.replace('Year','-year')
    os.system(f'cp -rnv ./{d}/* {args.out}/{args.subj}/{target_d}') # verbose, no overwriting
os.chdir(tmpdir)
os.system(f'rm -rf {args.subj}*')
os.chdir(args._in)
os.system(f'datalad drop {args.subj}')
# this removes the files but preserves the symbolic links, for further updates