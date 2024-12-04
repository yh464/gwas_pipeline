#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2024-11-14

A simplified script to conduct genetic correlation between groups of phenotypes

Requires following inputs: 
    GWAS summary statistics (scans directory for all files)
'''


def main(args):
    import os
    from fnmatch import fnmatch
    gcorrdir = f'{args.out}/gcorr'
    if not os.path.isdir(gcorrdir): os.system(f'mkdir -p {gcorrdir}')
    
    # array submitter
    from _utils import array_submitter
    submitter = array_submitter.array_submitter(
        name = f'gcorr_{args.p1[0]}_{args.p2[0]}',
        timeout = 10, mode = 'long')
    
    scripts_path = os.path.realpath(__file__)
    scripts_path = os.path.dirname(scripts_path)
    
    # scans directories to include sumstats 
    os.chdir(args._in)
    prefix_1 = []; pheno_1 = []
    prefix_2 = []; pheno_2 = []
    for p in args.p1:
        for x in os.listdir(p):
            if fnmatch(x,'*.sumstats'):
                prefix_1.append(x.replace('.sumstats','')); pheno_1.append(p)
    for p in args.p2:
        for x in os.listdir(p):
            if fnmatch(x,'*.sumstats'):
                prefix_2.append(x.replace('.sumstats','')); pheno_2.append(p)
    
    for i in range(len(prefix_1)):
      for j in range(len(prefix_2)):
        out_rg = f'{gcorrdir}/{pheno_1[i]}_{prefix_1[i]}.{pheno_2[j]}_{prefix_2[j]}.rg.log'
        out_rg1 = f'{gcorrdir}/{pheno_2[j]}_{prefix_2[j]}.{pheno_1[i]}_{prefix_1[i]}.rg.log'
        
        # check for duplicate files
        if os.path.isfile(out_rg) and os.path.isfile(out_rg1):
            os.remove(out_rg1)
        if not os.path.isfile(out_rg) and os.path.isfile(out_rg1):
            out_rg = out_rg1
        
        # QC out_rg file
        if os.path.isfile(out_rg):
          # if analysis is aborted
          if os.path.getsize(out_rg) < 2048:
            os.remove(out_rg)
        
        if os.path.isfile(out_rg):
          # if it runs NA 
          tmp = open(out_rg).read().splitlines()[-4].split()[2]
          try: float(tmp)
          except:
              # try constraining intercepts
            submitter.add('bash '+
              f'{scripts_path}/ldsc_master.sh ldsc.py --ref-ld-chr {args.ldsc}/baseline/'+
              f' --w-ld-chr {args.ldsc}/baseline/ '+
              f'--rg {args._in}/{pheno_1[i]}/{prefix_1[i]}.sumstats,'+
              f'{args._in}/{pheno_2[j]}/{prefix_2[j]}.sumstats '+
              f'--out {out_rg[:-4]} --no-intercept')
            continue
        
        if os.path.isfile(out_rg) and (not args.force): continue
        
        submitter.add('bash '+
          f'{scripts_path}/ldsc_master.sh ldsc.py --ref-ld-chr {args.ldsc}/baseline/'+
          f' --w-ld-chr {args.ldsc}/baseline/ '+
          f'--rg {args._in}/{pheno_1[i]}/{prefix_1[i]}.sumstats,'+
          f'{args._in}/{pheno_2[j]}/{prefix_2[j]}.sumstats '+
          f'--out {out_rg[:-4]}')
    
    submitter.submit()
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'This programme estimates genetic cross-correlation between two different sets of phenotypes')
    parser.add_argument('-p1', help = 'First group of phenotypes to correlate, usually IDP', nargs = '*')
    parser.add_argument('-p2', help = 'Second group of phenotypes to correlate, usually disorders', nargs = '*',
      default = ['global_structural','disorders','gradients'])
    parser.add_argument('-i','--in', dest = '_in', help = 'GWA file directory',
      default = '../gene_corr/ldsc_sumstats/')
    parser.add_argument('--ldsc', dest = 'ldsc', help = 'LDSC executable directory',
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/ldsc/') # intended to be absolute
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
      default = '../gene_corr/')
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
      default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','ldsc']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheng_%pheno_%maf.sumstats', __file__)
    proj.add_output(args.out+'/gcorr/%pheng_%pheno_%maf.%pheng_%pheno_%maf.rg.log', __file__)
    try: main(args)
    except: cmdhistory.errlog()