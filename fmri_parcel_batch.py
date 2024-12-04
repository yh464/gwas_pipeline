# -*- coding: utf-8 -*-
#!/usr/bin/env python3

'''
Authors:
Version 1: 2017-04-13
    Rafael Romero-Garcia rr480@cam.ac.uk
Version 2: 2024-11-17
    Yuankai He yh464@cam.ac.uk

A script to batch run parcellations on large fMRI datasets

Changelog:
    Add modifiable input options, including using the %sub wildcard
    Allows the selection of a subset of subjects using a list
    Modified to support fMRI data with multiple runs
'''

def main(args):
    import os
    from fnmatch import fnmatch
    from time import perf_counter as t
    
    # array submitter
    from _utils import array_submitter
    submitter = array_submitter.array_submitter(
        name = 'fmri_parcel', n_cpu = 3,
        timeout = 40, modules = ['fsl','matlab/r2021b','freesurfer'],
        debug = False
        )
    
    tic = t()
    
    # error log for subjects with missing data
    errlog = open('/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/logs/parcel_error.txt','w')
    
    # force parameter
    force = '-f' if args.force else ''
    
    os.chdir(args._dir)
    
    tmpid = 0
    subs = open(args.subj).read().splitlines() if os.path.isfile(args.subj) else os.listdir()
    for sub in subs:
        # substitute '%sub' in input data to simplify variable names
        surf = args._dir + '/' + args.surf.replace('%sub',sub)
        fmri = args._dir + '/' + args.fmri.replace('%sub',sub) # some subjects do not have fMRI at all
        if not os.path.isdir(fmri):
            print(f'{sub} no original fMRI data', file = errlog)
            # print(f'{sub} no original fMRI data')
            continue
        
        if not os.path.isdir(f'{fmri}/registration'): os.mkdir(f'{fmri}/registration')
        if not os.path.isdir(f'{fmri}/parcellations'): os.mkdir(f'{fmri}/parcellations')
        
        # check for FreeSurfer reconstruction
        if not os.path.isfile(f'{surf}/label/rh.aparc.annot'):
            print(f'{sub} FreeSurfer reconstruction missing', file = errlog)
            # print(f'{sub} FreeSurfer reconstruction missing')
            continue
        
        # check for preprocessed fMRI data
        flist = []
        for f in os.listdir(fmri):
            # if fnmatch(f, '*desc-preproc_bold.nii.gz'):
            if fnmatch(f, '*desc-smoothAROMAnonaggr_bold.nii.gz'):
                flist.append(f)
            # if fnmatch(f,'*bold.feat'):
            #     flist.append(f'{f}/filtered_func_data_clean.nii.gz')
        if len(flist) == 0:
            print(f'{sub} processed fMRI missing', file = errlog)
            # print(f'{sub} processed fMRI missing')
            continue
        flist.sort()
        
        # check for parcellation of structural space
        valid_parcs = []
        for parc in args.parc:
            if not os.path.isfile(f'{surf}/parcellation/{parc}.nii.gz'):
                # print(f'{sub} missing FreeSurfer parcellation for {parc}, skipping')
                print(f'{sub} missing FreeSurfer parcellation for {parc}, skipping', file = errlog)
                continue
            
            # make output folders
            valid_parcs.append(parc)
            parcout = f'{fmri}/parcellations/{parc}'
            if not os.path.isdir(parcout): os.mkdir(parcout)
        valid_parcs = ' '.join(valid_parcs)
        
        # parcellation and coregistration for each run
        for run in flist:
            # check progress
            tmp = run.split('_') # ..._run-01_desc-preproc_bold.nii.gz or ..._run-01_space-MNI..._desc..._bold.nii.gz
            runid = 'run-01'
            for x in tmp:
                if fnmatch(x, 'run*'):
                    runid = x; break
            
            done = True
            for parc in valid_parcs.split(' '):
                parcdir = f'{fmri}/parcellations/{parc}'
                out_ts = f'{parcdir}/{runid}_ts_sc2345.txt'
                out_mat = f'{parcdir}/{runid}_Connectivity_sc2345.txt'
                if not os.path.isfile(out_ts) or not os.path.isfile(out_mat):
                    done = False
                    
            if done and not args.force: continue
        
            submitter.add(f'python fmri_parcel.py {sub} {run} -d {args._dir} -c {args.code} '+
                          f'--fmri {args.fmri} -s {args.surf} -p {valid_parcs} {force}') 
        
        tmpid += 1
        if tmpid % 100 == 0:
            toc = t() - tic
            print(f'processed {tmpid} subjects, time = {toc:.3f} seconds')
    submitter.submit()
    
if __name__ == '__main__':
    import argparse
    import os
    parser = argparse.ArgumentParser(
        description = 'This programme runs parcellation for all subjects across a dataset')
    parser.add_argument('-i','--subj', dest = 'subj', help = 'List of subjects - if not a file, then automatically scans dir',
        default = '../params/subjlist_abcd_eur.txt')
    parser.add_argument('-d','--dir', dest = '_dir', help = 'Directory containing imaging data',
        default = '/rds/project/rb643/rds-rb643-abcd/Data_Imaging') # intended to be absolute
    parser.add_argument('-c','--code', dest = 'code', help = 'Directory of matlab code for parcellation',
        default = os.path.realpath('../fmri'))
    parser.add_argument('--fmri', dest = 'fmri', 
        help = 'Path from data folder to fMRI data, %sub = placeholder for subject',
        default = '%sub/ses-baseline-year1/func')
    parser.add_argument('-s','--surf', dest = 'surf', help = 'Path from data folder to freesurfer output',
        default = '%sub/ses-baseline-year1/surfaces/%sub')
    parser.add_argument('-p','--parc', dest = 'parc', nargs = '*', help = 'Parcellations',
        #default = ['aparc_seq','500_sym.aparc_seq','HCP.fsaverage.aparc_seq','sjh_seq','economo_seq']
        default = ['HCP.fsaverage.aparc_seq']
        )
    parser.add_argument('-f','--force',dest = 'force', help = 'force overwrite',
        default = False, action = 'store_true')
    args = parser.parse_args()
    
    from _utils import cmdhistory
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()