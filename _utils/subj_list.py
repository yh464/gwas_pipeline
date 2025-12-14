'''
This script screens subjects with a valid imaging profile
'''

from .logger import logger
log = logger()

def main(args):
    import os
    
    # fail-safe
    args.out = os.path.realpath(args.out)
    if not os.path.isdir(args.out): os.mkdir(args.out)
    
    # progress check
    fout = f'{args.out}/{args.prefix}.txt'
    errlog = f'{args.out}/{args.prefix}_not_found.txt'
    if os.path.isfile(fout) and not args.force: 
        log.log('subj list already generated')
        return
    
    # count subjs with imaging profiles and w/o
    found = 0
    not_found = 0
    fout = open(fout,'w')
    errlog = open(errlog,'w')
    base = args.target.split('%subj')[0]
    for subj in os.listdir(base):
        target = args.target.replace('%subj',subj) # target file path
        if os.path.isfile(target):
            found += 1
            log.log(subj.replace('UKB',''), file = fout)
        else:
            not_found += 1
            log.log(subj.replace('UKB',''), file = errlog)
    
    log.log(f'Total {found + not_found} subjects')
    log.log(f'Found imaging profiles for {found} subjects')
    log.log(f'No imaging profile for {not_found} subjects')
    
    return

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='This programme finds subjects with '+
        'a valid imaging profile')
    parser.add_argument('-t','--target',dest = 'target', help =
        'Target file to screen',
        default = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Imaging/'+
        '%subj/func/fMRI/parcellations/HCP.fsaverage.aparc_seq/Connectivity_sc2345.txt')
    parser.add_argument('-o','--out', dest = 'out', help = 'output subj list dir',
        default = '../params')
    parser.add_argument('-p','--prefix', dest = 'prefix', required = True,
        help = 'output file prefix')
    parser.add_argument('-f','--force', dest = 'force', action = 'store_true',
        default = False, help = 'force overwrite')
    args = parser.parse_args()
    
    from _utils import cmdhistory
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()