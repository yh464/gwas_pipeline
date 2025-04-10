#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2023-07-13
Version 2.0: 2024-11-12
Version 2.1: 2024-11-29

Conducts genetic correlation for groups of phenotypes,
including correlations within and betwwen groups of phenotypes

Preceding workflow:
    heri_batch.py
Requires following inputs:
    MUNGED GWAS summary statistics (scans directory for all files)
Changelog:
    moved all heritability/munging elements to heri_batch.py
    Added the feature to check for duplicate files (a with b v b with a)
'''

def main(args):
    import os
    from fnmatch import fnmatch
    if not os.path.isdir(args.out): os.system(f'mkdir -p {args.out}')
    
    # array submitter
    from _utils import array_submitter
    submitter = array_submitter.array_submitter(
        name = 'gcorr_'+'_'.join(args.pheno),
        timeout = 10, mode = 'long',
        debug = False
        )
    
    scripts_path = os.path.realpath(__file__)
    scripts_path = os.path.dirname(scripts_path)
    prefix_list = []
    pheno_list = []
    
    for x in args.pheno:
      for y in os.listdir(f'{args._in}/{x}'):
        if fnmatch(y, '*.sumstats'):
          prefix = y.replace('.sumstats','')
          prefix_list.append(prefix)
          pheno_list.append(x)
          
    os.chdir(args._in)
    
    for i in range(len(prefix_list)):
      for j in range(i):
        g1 = pheno_list[i]; p1 = prefix_list[i]
        g2 = pheno_list[j]; p2 = prefix_list[j]
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
                f'--rg {args._in}/{g1}/{p1}.sumstats,'+
                f'{args._in}/{g2}/{p2}.sumstats '+
                f'--out {out_rg[:-4]} --no-intercept')
              continue
        
        # if out_rg survives QC we continue without doing anything
        if os.path.isfile(out_rg) and (not args.force): continue
        
        # if no out_rg file is present
        submitter.add('bash '+
          f'{scripts_path}/ldsc_master.sh ldsc.py --ref-ld-chr {args.ldsc}/baseline/'+
          f' --w-ld-chr {args.ldsc}/baseline/ '+
          f'--rg {args._in}/{g1}/{p1}.sumstats,'+
          f'{args._in}/{g2}/{p2}.sumstats '+
          f'--out {out_rg[:-4]}')
    
    submitter.submit()
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 
      'This programme creates genetic correlation matrices for global phenotypes')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*')
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
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_output(args._in+'/ldsc_sumstats/%pheng_%pheno_%maf.sumstats', __file__)
    proj.add_output(args.out+'/rglog/%pheng.%pheng/%pheng_%pheno.%pheng_%pheno.rg.log', __file__)
    try: main(args)
    except: cmdhistory.errlog()