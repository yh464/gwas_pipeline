import os
wd = '/rds/project/rds-CeXlNYOYMxw/userdata/yh464/ABCD_fMRIprep/fmriprep/'
os.chdir(wd)

from _utils import array_submitter
submitter = array_submitter.array_submitter(
    name = 'abcd_download', n_cpu = 1,
    timeout = 10, mode = 'long', arraysize = 100,
    env = 'datalad', wd = wd)

dirlist = open('../subj_list.txt').read().splitlines()
for i in range(len(dirlist)):
   submitter.add(f'datalad get {dirlist[i]}/ses*/func/*AROMA* {dirlist[i]}/ses*/func/*desc-confounds_timeseries.tsv {dirlist[i]}/ses*/func/*desc-preproc_bold.nii.gz '+
                 f'{dirlist[i]}/ses*/anat/*desc-preproc_T1w.nii.gz {dirlist[i]}/ses*/anat/*surf.gii')
submitter.submit()
