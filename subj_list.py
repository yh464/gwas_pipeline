'''
This script screens subjects with a valid imaging profile
'''

def main(args):
    import os
    
    # fail-safe
    outdir = os.path.dirname(args.out)
    if not os.path.isdir(outdir): os.mkdir(outdir)
    
    # progress check
    fout = f'{args.out}.txt'
    errlog = args.out.replace('.txt','_not_found.txt')
    if os.path.isfile(fout) and not args.force: 
        print('subj list already generated')
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
            print(subj.replace('UKB',''), file = fout)
        else:
            not_found += 1
            print(subj.replace('UKB',''), file = errlog)
    
    print(f'Total {found + not_found} subjects')
    print(f'Found imaging profiles for {found} subjects')
    print(f'No imaging profile for {not_found} subjects')
    return

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='This programme finds subjects with '+
        'a valid imaging profile')
    parser.add_argument('-t','--target',dest = 'target', help =
        'Target file to screen',
        default = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Imaging/'+
        '%subj/func/fMRI/parcellations/HCP.fsaverage.aparc_seq/Connectivity_sc2345.txt')
    parser.add_argument('-o','--out', dest = 'out', help = 'output prefix', required = True)
    parser.add_argument('-f','--force', dest = 'force', action = 'store_true',
        default = False, help = 'force overwrite')
    args = parser.parse_args()
    
    from ._utils import cmdhistory
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()