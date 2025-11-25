#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-01-20

Utility script to extract input SNPs from GWAS files

Requires following inputs: 
    GWAS summary statistics (single file)
'''

def search_file(file, patterns):
    import os
    import numpy as np
    import pandas as pd
    import subprocess
    
    patterns = os.path.realpath(patterns)
    hdr = open(file).readline().replace('\n','').split()
    cmd = ['/bin/fgrep','-wf',patterns,file]
    search = subprocess.Popen(cmd, stdout = subprocess.PIPE)
    
    try: df = pd.read_table(search.stdout, header = None)
    except: print(f'None of the SNPs found in {file}'); return None
    
    df.columns = hdr
    df['Phenotype'] = os.path.basename(file).replace('.fastGWA','')
    if 'OR' in df.columns: df['BETA'] = np.log(df['OR'])
    if 'N' not in df.columns: df['N'] = df.N_CAS + df.N_CON
    for col in ['BETA', 'SE', 'AF1']:
        if col not in df.columns: df[col] = np.nan
    df = df.loc[:,['Phenotype','SNP','BETA','SE','P','N','A1','A2','AF1','CHR','POS']]
    return df
    
def search_snp(x, tmpdir, args): 
    # must define in main environment for multiprocessing
    import os
    from fnmatch import fnmatch
    import numpy as np
    import pandas as pd
    from multiprocessing import Pool
    from functools import partial
    print(x)
    
    if os.path.isfile(args.snp[0]) and len(args.snp) == 1:
        patterns = args.snp[0]
        snps = open(patterns).read().splitlines()
    else: patterns = f'{tmpdir}/snp_list.txt'; snps = args.snp
    # search for required SNP using grep
    flist = []
    for y in os.listdir(f'{args._in}/{x}'):
        if fnmatch(y,'*.fastGWA') and not fnmatch(y, '*_X.fastGWA'): 
            flist.append(f'{args._in}/{x}/{y}')
    if len(flist) == 0: print(f'NO GWAS FILE FOR {x}'); return None
    
    cache = f'{tmpdir}/sigsnp_{snps[0]}_{snps[-1]}_{x}.txt'
    if os.path.isfile(cache) and not args.force:
        df = pd.read_table(cache)
        return df
    
    pool = Pool(min((len(flist)),10))
    search = pool.map(partial(search_file, patterns=patterns), flist, 
                      chunksize = int(np.ceil(len(flist)/10)))
    temp = []
    for y in search:
        if type(y) != type(None): temp.append(y)
    df = pd.concat(temp)
    df.insert(loc = 0, column = 'Group', value = x)
    df.insert(loc = 2, column = 'q', value = df.P * 1e+6)
    df.insert(loc = 4, column = 'Z', value = df.BETA/df.SE)
    df.to_csv(cache, sep = '\t', index = False)
    return df

def main(args):
    import os
    import pandas as pd
    from ._utils.path import normaliser
    
    norm = normaliser()
    all_files = []
    tmpdir = os.path.realpath('../temp/single_snp')
    if not os.path.isdir(tmpdir): os.system(f'mkdir -p {tmpdir}')
    
    if os.path.isfile(args.snp[0]) and len(args.snp) == 1:
        patterns = args.snp[0]
    else:
        patterns = f'{tmpdir}/snp_list.txt'
        with open(patterns,'w') as tmpfile:
            for snp in args.snp: print(snp, file = tmpfile)
            tmpfile.close()
    
    temp = [search_snp(x, tmpdir, args) for x in args.pheno]
    all_files = []
    for x in temp: 
        if type(x) != type(None): all_files.append(x)
    
    all_files = pd.concat(all_files).sort_values(by = ['Group','Phenotype','SNP'])
    if args.out != None: 
        all_files.to_csv(args.out, sep = '\t', index = False)
    norm.normalise(all_files).to_clipboard(index = False)
    # all_files_wide = all_files.pivot_table(values = 'p', index = 'SNP', columns = 'Phenotype')
    # all_files_wide.insert(loc = 0, column = 'min_pval', value = all_files_wide.min(axis = 1))
    # print(norm.normalise(all_files_wide))
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
        'This script extract selected SNPs from GWAS sumstats')
    parser.add_argument('snp',help = 'list of SNPs to extract, txt file or as an argument',
        nargs = '*')
    parser.add_argument('-p', dest = 'pheno', help = 'list of phenotype groups to extract SNPs',
        nargs = '*', default = ['global_graph','global_asym'])
    parser.add_argument('-i','--in', dest = '_in', 
        help = 'Directory containing all GWA summary statistics',
        default = '../gwa/')
    parser.add_argument('-o','--out', dest = 'out', help = 'Output list file')
    parser.add_argument('-f','--force', dest = 'force', action = 'store_true',
        default = False, help = 'force overwrite')
    args = parser.parse_args()
    import os
    args._in = os.path.realpath(args._in)
    if args.out != None: args.out = os.path.realpath(args.out)
    args.pheno.sort()
    
    from ._utils import cmdhistory, logger
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()