#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-02-18

A utility to rename files across the project
'''

protected_directories = [
    'path',
    'cmdhistory',
    'toolbox',
    'scripts',
    'params',
    'fmri',
    'temp',
    'archive-2023-rsfmri-gwas'
    ]

def re(d, before, after, protect):
    import os
    from fnmatch import fnmatch
    if os.path.basename(d) in protected_directories: return
    print(d)
    if len(os.listdir(d)) > 50000:
        print(f'Skipping {d} due to large number of files')
        return
    for x in os.listdir(d):
        if os.path.isdir(f'{d}/{x}'):
            re(f'{d}/{x}', before, after, protect)
        if any([fnmatch(x,f'*{p}*') for p in protect]): continue
        y = x
        for b, a in zip(before, after):
            y = y.replace(b,a)
        os.rename(f'{d}/{x}', f'{d}/{y}')

def main(args):
    re(args._dir, args.before, args.after, args.protect)
    return

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Computes correlational plots between IDPs and PGS (batch)')
    parser.add_argument('before', nargs = '*', help = 'File names to rename')
    parser.add_argument('-a','--after',nargs = '*', help = 'File names after renaming, corresponding to "before"', 
                        required = True)
    parser.add_argument('-p','--protect', nargs = '*', help = 'Protected file names that will not be changed',
                        default = [])
    parser.add_argument('-d','--dir', dest = '_dir', help = 'Project root directory',
        default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/')
    args=parser.parse_args()
    import os
    args._dir = os.path.realpath(args._dir)

    import cmdhistory, logger
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()