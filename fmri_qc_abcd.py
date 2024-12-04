#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1.0: 2024-11-21
Version 1.1: 2024-11-26

Quality controls fMRI data in the ABCD dataset

Preceding workflow:
    fmri_parcel_batch.py
Requires following inputs: 
    fMRI data (to be scanned)
    desc confounds time series
    connectivity matrix
'''
def parse_aroma_cov(conf, sub, run):
    import pandas as pd
    df = pd.read_table(conf)
    nframes = df.shape[0]
    fd = df.loc[1:,'framewise_displacement'] # no FD for first time point
    fdmax = fd.max()
    fd = fd.mean()
    rmsd = df.loc[1:,'rmsd'] # no RMSD for first time point
    rmsdmax = rmsd.max()
    rmsd = rmsd.mean()
    outdf = pd.DataFrame(dict(
        sub = sub,
        run = run,
        length = nframes,
        fd = fd, fdmax = fdmax,
        rmsd = rmsd, rmsdmax = [rmsdmax] # list format is needed so length is one
        ))
    return outdf

def parse_connectivity_mat(conn, sub, run):
    import numpy as np
    import pandas as pd
    # from scipy.spatial.distance import squareform
    con = np.loadtxt(conn, delimiter = ',')
    # con = 0.5 * con + 0.5 * con.T
    tmp = []
    for i in range(con.shape[0]):
        for j in range(i):
            tmp.append(con[i,j])
    con = np.array(tmp)
    # con = squareform(con)
    m = con.mean()
    l = con.min()
    h = con.max()
    nnodes = con.shape[0]
    outdf = pd.DataFrame(dict(sub = [sub], run = run, mean_conn = m, min_conn = l, max_conn = h, nnodes = nnodes))
    return outdf

# aggregates connectivity matrices for the same subject
def avg_connectivity_mat(mats, weights, out):
    import numpy as np
    l = []
    for m in mats:
        tmp = np.loadtxt(m, delimiter = ',')
        tmp = (tmp + tmp.T)/2
        for i in range(tmp.shape[0]): tmp[i,i] = 0
        tmp = np.arctanh(tmp)
        l.append(tmp)
    try:
        l = np.average(l, axis = 0, weights = weights)
        l = np.tanh(l)
        np.savetxt(out, l, delimiter = ',')
    except: print(mats.tolist())

def mad_threshold(x, d = 5):
    import numpy as np
    x = np.array(x)
    m = np.median(x)
    ad = np.abs(x - m)
    mad = np.median(ad)
    return [m - d * mad, m + d * mad]
    
def main(args):
    import os
    from fnmatch import fnmatch
    import pandas as pd
    from time import perf_counter as t
    import numpy as np
    
    tic = t()
    if not os.path.isfile(f'{args.out}/abcd_fmri_qc.txt') or args.force:
        os.chdir(args._dir)
        subs = open(args.subj).read().splitlines() if os.path.isfile(args.subj) else os.listdir()
        all_conf = []
        all_conn = []
        nsub = 0
        for sub in subs:
            fmri = args._dir + '/' + args.fmri.replace('%sub',sub) # some subjects do not have fMRI at all
            if not os.path.isdir(fmri):
                print(f'{sub} no original fMRI data')
                continue
            
            # check for preprocessed fMRI data
            for f in os.listdir(fmri):
                # if fnmatch(f, '*desc-preproc_bold.nii.gz'):
                if fnmatch(f, '*desc-smoothAROMAnonaggr_bold.nii.gz'):
                    confounds = fmri + '/' + f.replace(
                        '_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii.gz',
                        '_desc-confounds_timeseries.tsv')
                    if not os.path.isfile(confounds): continue
                    
                    # collate connectivity matrices
                    for tmp in f.split('_'):
                        if fnmatch(tmp, 'run-*'): runid = tmp; break
                    conn = f'{fmri}/parcellations/{args.parc}/{runid}_Connectivity_sc2345.txt'
                    if not os.path.isfile(conn): continue
                    
                    # collate statistics
                    stats = parse_aroma_cov(confounds, sub, conn)
                    all_conf.append(stats)
                    
                    con = parse_connectivity_mat(conn, sub, conn)
                    all_conn.append(con)                
            nsub += 1
            if nsub % 100 == 0:
                toc = t() - tic
                print(f'{nsub} subjects loaded, time = {toc:.3f} seconds')
        all_conf = pd.concat(all_conf)
        all_conn = pd.concat(all_conn)
        
        all_runs = pd.merge(all_conf, all_conn, on = ['sub', 'run'])
        all_runs.to_csv(f'{args.out}/abcd_fmri_qc.txt', sep = '\t', index = False)
    
    all_runs = pd.read_table(f'{args.out}/abcd_fmri_qc.txt')
    # drop runs with insufficient length
    all_runs = all_runs.loc[all_runs.length >= 100,:]
    
    # collate 5MAD thresholds for all stats
    fd_thr = mad_threshold(all_runs.fd)[1] # only +ve threshold
    fdmax_thr = mad_threshold(all_runs.fdmax)[1]
    rmsd_thr = mad_threshold(all_runs.rmsd)[1]
    rmsdmax_thr = mad_threshold(all_runs.rmsdmax)[1]
    
    _filter = (all_runs.fd < fd_thr) * (all_runs.fdmax < fdmax_thr) * \
        (all_runs.rmsd < rmsd_thr) * (all_runs.rmsdmax < rmsdmax_thr) * \
        (all_runs.nnodes == np.median(all_runs.nnodes))
    _filter = _filter.astype('?')
    
    all_runs = all_runs.loc[_filter, :]
    subs = all_runs['sub'].unique()
    print(f'{all_runs.shape[0]} runs remaining from {len(subs)} subjects')
    print(f'mean connecivity {all_runs.mean_conn.mean():.4f}')
    print(f'min connectivity mean {all_runs.min_conn.mean():.4f}')
    all_runs.to_csv(f'{args.out}/abcd_fmri_qced.txt', sep = '\t', index = False)
    
    # export counts of runs retained after QC as a covariate
    count = []
    for s in subs:
        count.append(all_runs.loc[all_runs['sub'] == s, :].shape[0])
    pd.DataFrame(dict(sub = subs, n_runs = count)).to_csv(
        f'{args.out}/abcd_fmri_runcount.txt', sep = '\t', index = False)
    
    # aggregate connectivity matrices for individual subjects
    nsub = 0
    for s in subs:
        fmri = args._dir + '/' + args.fmri.replace('%sub',s) # some subjects do not have fMRI at all
        out = f'{fmri}/parcellations/{args.parc}/Connectivity_sc2345.txt'
        # if not os.path.isfile(out) or args.force:
        if True:
            tmp = all_runs.loc[all_runs['sub'] == s, ['run', 'length']]
            avg_connectivity_mat(tmp['run'], tmp['length'],out)
        nsub += 1
        if nsub % 100 == 0:
            toc = t() - tic
            print(f'{nsub} subjects completed aggregation of matrices, time = {toc:.3f} seconds')
    return

if __name__ == '__main__':
    import argparse
    import os
    parser = argparse.ArgumentParser(
        description = 'This programme runs parcellation for all subjects across a dataset')
    parser.add_argument('-i','--subj', dest = 'subj', help = 'List of subjects - if not a file, then automatically scans dir',
        default = '../params/subjlist_abcd_eur.txt')
    parser.add_argument('-d','--dir', dest = '_dir', help = 'Root directory containing imaging data',
        default = '/rds/project/rb643/rds-rb643-abcd/Data_Imaging') # intended to be absolute
    parser.add_argument('--fmri', dest = 'fmri', 
        help = 'Path from data folder to fMRI data, %sub = placeholder for subject',
        default = '%sub/ses-baseline-year1/func')
    parser.add_argument('-p','--parc', dest = 'parc', help = 'Parcellation to QC',
        default = 'HCP.fsaverage.aparc_seq')
    parser.add_argument('-o', '--out', dest = 'out', help = 'output QC file',
        default = '../params')
    parser.add_argument('-f','--force',dest = 'force', help = 'force overwrite',
        default = False, action = 'store_true')
    args = parser.parse_args()
    
    args.subj = os.path.realpath(args.subj)
    args.out = os.path.realpath(args.out)
    
    from _utils import cmdhistory
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()