#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2024-11-14
Version 2: 2025-04-09

A simplified script to conduct genetic correlation between groups of phenotypes

Requires following inputs: 
    GWAS summary statistics (scans directory for all files)
'''


def main(args):
    import os
    from fnmatch import fnmatch
    if not os.path.isdir(args.out): os.system(f'mkdir -p {args.out}')
    
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
        g1 = pheno_1[i]; p1 = prefix_1[i]
        g2 = pheno_2[j]; p2 = prefix_2[j]
        if g2 < g1: g1,g2 = (g2,g1); p1,p2 = (p2,p1)
        
        out_rg = f'{args.out}/{g1}.{g2}/{g1}_{p1}.{g2}_{p2}.rg.log'
        if not os.path.isdir(f'{args.out}/{g1}.{g2}'): 
            os.mkdir(f'{args.out}/{g1}.{g2}')
        
        # check for duplicate files
        for out_rg1 in [f'{args.out}/{g1}.{g2}/{g1}_{p1}.{g2}_{p2}.rg.log',
                        f'{args.out}/{g1}.{g2}/{g2}_{p2}.{g1}_{p1}.rg.log',
                        f'{args.out}/{g2}.{g1}/{g1}_{p1}.{g2}_{p2}.rg.log',
                        f'{args.out}/{g2}.{g1}/{g2}_{p2}.{g1}_{p1}.rg.log']:
            if os.path.isfile(out_rg) and os.path.isfile(out_rg1) and out_rg1 != out_rg:
                os.remove(out_rg1)
            if not os.path.isfile(out_rg) and os.path.isfile(out_rg1):
                os.rename(out_rg1, out_rg)
        
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
      default = '../gcorr/ldsc_sumstats/')
    parser.add_argument('--ldsc', dest = 'ldsc', help = 'LDSC executable directory',
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/ldsc/') # intended to be absolute
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
      default = '../gcorr/rglog/')
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
      default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','ldsc']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
    from _utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheng_%pheno_%maf.sumstats', __file__)
    proj.add_output(args.out+'/%pheng_%pheno_%maf.%pheng_%pheno_%maf.rg.log', __file__)
    try: main(args)
    except: cmdhistory.errlog()