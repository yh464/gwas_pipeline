#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
2024-10-25

Preprocesses fMRI images as per the UKBiobank pipeline
'''

if __name__ == '__main__':
    import argparse
    import os
    parser = argparse.ArgumentParser(
        description = 'This programme runs parcellation for all subjects across a dataset')
    
    parser.add_argument('sub', help = 'Subject ID')
    parser.add_argument('img', help = 'Name of fMRI image file, under fmri directory below')
    
    parser.add_argument('-d','--dir', dest = '_dir', help = 'Directory containing imaging data',
        default = '/rds/project/rb643/rds-rb643-abcd/Data_Imaging') # intended to be absolute
    parser.add_argument('-c','--code', dest = 'code', help = 'Directory of design file for FEAT analysis',
        default = '../fmri/design_abcd.fsf')
    
    parser.add_argument('--fmri', dest = 'fmri', 
        help = 'Path from data folder to fMRI data, %sub = placeholder for subject',
        default = '%sub/ses-baseline-year1/func')
    parser.add_argument('-s','--surf', dest = 'surf', help = 'Path from data folder to freesurfer output',
        default = '%sub/ses-baseline-year1/surfaces/%sub')
    
    parser.add_argument('-f','--force',dest = 'force', help = 'force overwrite',
        default = False, action = 'store_true')
    args = parser.parse_args()
    args.code = os.path.realpath(args.code)
    
def main(args):
    import os
    from fnmatch import fnmatch
    
    # directory management
    surf = f'{args._dir}/{args.surf}'.replace('%sub',args.sub)
    fmri = f'{args._dir}/{args.fmri}'.replace('%sub',args.sub)
    
    # image prefix
    prefix = args.img.replace('.nii','').replace('.gz','')
    prefix = os.path.basename(prefix)
    tmp = prefix.split('_')
    for tmp1 in tmp:
        if fnmatch(tmp1,'run*'): runid = tmp1
    del tmp, tmp1
    
    # convert structural image
    t1orig = f'{surf}/mri/T1_orig.nii.gz'
    if not os.path.isfile(t1orig) or args.force:
        os.system(f'mri_convert {surf}/mri/orig.mgz {t1orig}')
    
    # extract single band reference image
    sbref = f'{fmri}/rfMRI.SBREF'
    if not os.path.isdir(sbref): os.mkdir(sbref)
    sbref_img = f'{sbref}/{prefix}_SBREF.nii'
    if not os.path.isfile(sbref_img):
        os.system(f'3dcalc -a {fmri}/{args.img} -expr a -prefix {sbref_img}')
    
    # edit FEAT design file to match the directories, then execute FEAT using the UKB pipeline except B0 unwarping
    featdir = f'{fmri}/{prefix}.feat'
    if not os.path.isfile(f'{featdir}/filtered_func_data.nii.gz') or args.force:
        os.system(f'rm -rf {featdir}') # this is to ensure that FEAT output is consistent
        design = f'{fmri}/design_{runid}.fsf'
        os.system(f'cp -f {args.code} {design}')
        f = open(design,'a')
        print(f'set feat_files(1) \"{fmri}/{prefix}\"', file = f)
        print(f'set alt_ex_func(1) \"{sbref_img}\"', file = f)
        print(f'set highres_files(1) \"{t1orig}\"', file = f)
        print(f'set fmri(outputdir) \"{featdir}\"', file = f)
        f.close()
        os.system(f'feat {design}')
    
    # fix to obtain ICA filtered data
    if not os.path.isfile(f'{featdir}/filtered_func_data_clean.nii.gz') or args.force:
        os.system(f'fix {featdir} UKBiobank 20 -m')

main(args)