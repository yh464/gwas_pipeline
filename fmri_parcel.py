# -*- coding: utf-8 -*-
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

if __name__ == '__main__':
    import argparse
    import os
    parser = argparse.ArgumentParser(
        description = 'This programme runs parcellation for all subjects across a dataset')
    
    parser.add_argument('sub', help = 'Subject ID')
    parser.add_argument('img', help = 'Name of fMRI image file, NO directory')
    
    parser.add_argument('-d','--dir', dest = '_dir', help = 'Directory containing imaging data',
        default = '/rds/project/rb643/rds-rb643-abcd/Data_Imaging') # intended to be absolute
    parser.add_argument('-c','--code', dest = 'code', help = 'Directory of matlab code for parcellation',
        default = os.path.realpath('../fmri'))
    
    parser.add_argument('--fmri', dest = 'fmri', 
        help = 'Path from data folder to fMRI data, %sub = placeholder for subject',
        default = '%sub/ses-baseline-year1/func')
    parser.add_argument('-s','--surf', dest = 'surf', help = 'Path from data folder to freesurfer output',
        default = '%sub/ses-baseline-year1/surfaces/%sub')
    parser.add_argument('-p','--parc', dest = 'parc', nargs = '*', help = 'Parcellations to run',
        required = True)
    
    parser.add_argument('-f','--force',dest = 'force', help = 'force overwrite',
        default = False, action = 'store_true')
    args = parser.parse_args()
    
def main(args):
    # file names
    from fnmatch import fnmatch
    tmp = args.img.split('_') # ..._run-01_desc-preproc_bold.nii.gz or ..._run-01_space-MNI..._desc..._bold.nii.gz
    runid = 'run-01'
    for x in tmp:
        if fnmatch(x, 'run*'):
            runid = x; break
    surf = f'{args._dir}/{args.surf}'.replace('%sub',args.sub)
    fmri = f'{args._dir}/{args.fmri}'.replace('%sub',args.sub)
    reg = f'{fmri}/registration' # registration directory
    if not os.path.isdir(reg): os.mkdir(reg)
    
    # check progress
    skip = True
    for parc in args.parc:
        parcdir = f'{fmri}/parcellations/{parc}'
        if not os.path.isdir(parcdir): os.mkdir(parcdir)
        
        ## file names specific to parcellation
        out_ts = f'{parcdir}/{runid}_ts_sc2345.txt'
        out_mat = f'{parcdir}/{runid}_Connectivity_sc2345.txt'
        if not os.path.isfile(out_ts) or not os.path.isfile(out_mat): skip = False
    if skip and not args.force: 
        # clear up
        img_std = f'{reg}/fMRI_{runid}_std.nii.gz'
        img_t1space = f'{reg}/fMRI_{runid}_t1space.nii.gz'
        os.system(f'rm -f {img_std} {img_t1space}')
        return
    
    img = f'{fmri}/{args.img}'
    t12std_xfm = f'{reg}/T1_to_std.txt'
    fmri2t1_xfm = f'{reg}/fMRI_{runid}_to_T1.txt'
    t12fmri_xfm = f'{reg}/T1_to_fMRI_{runid}.txt'
    t1_std = f'{reg}/T1_std.nii.gz'
    img_std = f'{reg}/fMRI_{runid}_std.nii.gz'
    img_t1space = f'{reg}/fMRI_{runid}_t1space.nii.gz'
    
    print('='*20)
    print(f'Running co-registration and parcellation for subject {args.sub}, {runid}')
    
    # quality ctrl to remove broken files
    # import numpy as np
    # for tmp in [t12fmri_xfm, t12std_xfm]:
    #     if os.path.isfile(tmp):
    #         try:
    #             test = np.loadtxt(tmp, delimiter = '\s+')
    #             if test.shape != (4,4): os.remove(tmp)
    #             del test
    #         except:
    #             os.remove(tmp)
    # del tmp
    
    # first do co-registration
    print('Conducting co-registration')
    ## orientation
    if not os.path.isfile(f'{surf}/mri/T1_from_mask.nii.gz') or args.force:
        os.system(f'mri_convert {surf}/mri/brainmask.mgz {surf}/mri/T1_from_mask.nii.gz')
    
    if not os.path.isfile(t1_std) or not os.path.isfile(t12std_xfm) or args.force:
        os.system(f'fslreorient2std {surf}/mri/T1_from_mask.nii.gz {t1_std}')
        os.system(f'fslreorient2std {surf}/mri/T1_from_mask.nii.gz > {t12std_xfm}')
    
    if os.path.isfile(img_std) and args.force:
        os.remove(img_std)
    if not os.path.isfile(img_std):
        os.system(f'fslreorient2std {img} {img_std}')
        
    ## co-registration with FLIRT
    if not os.path.isfile(img_t1space) or args.force:
        os.system(f'flirt -in {img_std} -noresample -ref {t1_std} -out {img_t1space} -omat {fmri2t1_xfm}');
        os.system(f'convert_xfm {fmri2t1_xfm} -inverse -omat {t12fmri_xfm}')
    
    # then do parcellation
    for parc in args.parc:
        print(f'Conducting parcellation for {parc}')
        parcdir = f'{fmri}/parcellations/{parc}'
        if not os.path.isdir(parcdir): os.mkdir(parcdir)
        
        ## file names specific to parcellation
        parcimg = f'{surf}/parcellation/{parc}.nii.gz'
        parc_fmrispace = f'{parcdir}/{parc}_fMRI_{runid}_space.nii.gz'
        parc_std = f'{parcdir}/{parc}_std.nii.gz'
        out_ts = f'{parcdir}/{runid}_ts_sc2345.txt'
        out_mat = f'{parcdir}/{runid}_Connectivity_sc2345.txt'
        
        ## QC
        if os.path.isfile(parc_fmrispace):
            if os.path.getsize(parc_fmrispace) < 5000: os.remove(parc_fmrispace)
        
        ## co-register parcellation with 
        if not os.path.isfile(parc_fmrispace) or args.force:
            os.system(f'flirt -in {parcimg} -ref {t1_std} -applyxfm -init {t12std_xfm} -out {parc_std}')
            os.system(f'flirt -in {parc_std} -ref {img_std} -applyxfm -init {t12fmri_xfm} '+
                      f'-interp nearestneighbour -out {parc_fmrispace}')
        
        ## extract time series
        if not os.path.isfile(out_ts) or not os.path.isfile(out_mat) or args.force:
            os.system(f'matlab -nodisplay -nosplash -sd {args.code} -batch "fMRI_extraction(\'{img_std}\','+
                      f'\'{parc_fmrispace}\',\'{parcdir}/\',\'{runid}\')"')
            # -sd: starting directory, should contain both fMRI_extraction and wavelets files
    
    # clear up
    os.system(f'rm -f {img_std} {img_t1space}')
    
main(args)