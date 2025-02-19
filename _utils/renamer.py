#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-01-10

A utility to rename files across the project
'''

protected_directories = [
    'toolbox',
    'scripts',
    'fmri',
    'temp'
    ]

def re(d, before, after):
    import os
    if os.path.basename(d) in protected_directories: return
    print(d)
    if len(os.listdir(d)) > 50000:
        print(f'Skipping {d} due to large number of files')
        return
    for x in os.listdir(d):
        if os.path.isdir(f'{d}/{x}'):
            re(f'{d}/{x}', before, after)
        for b, a in zip(before, after):
            y = x.replace(b,a)
            os.rename(f'{d}/{x}', f'{d}/{y}')

def main(args):
    re(args._dir, args.before, args.after)
    return

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Computes correlational plots between IDPs and PGS (batch)')
    parser.add_argument('before', nargs = '*', help = 'File names to rename')
    parser.add_argument('-a','--after',nargs = '*', help = 'File names after renaming, corresponding to "before"')
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