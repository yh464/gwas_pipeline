# This script merges the ABCD dataset

subjlist = open('/rds/project/rds-CeXlNYOYMxw/userdata/yh464/ABCD_fMRIprep/subj_list.txt').read().splitlines()

from _utils import array_submitter
submitter = array_submitter.array_submitter(
    name = 'abcd_merge', n_cpu = 1,
    timeout = 10, mode = 'long', arraysize = 100,
    env = 'datalad')

for i in subjlist:
   submitter.add(f'python abcd_merge.py {i}')

submitter.submit()
