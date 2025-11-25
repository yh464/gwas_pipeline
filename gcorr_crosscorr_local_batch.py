#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-02-11

A simplified script to conduct genetic correlation between regional phenotypes 
and external correlates

Requires following inputs: 
    GWAS summary statistics (scans directory for all files)
'''


def main(args):
    import os
    from fnmatch import fnmatch
    if not os.path.isdir(args.out): os.system(f'mkdir -p {args.out}')
    
    # array submitter
    from _utils.slurm import array_submitter
    submitter = array_submitter(
        name = f'gcorr_{args.p1[0]}_{args.p2[0]}',
        timeout = 10, mode = 'long',
        debug = True
        )
    
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
    
    for g1, p1 in zip(pheno_1, prefix_1):
      if not os.path.isdir(f'{args.out}/{g1}'): os.mkdir(f'{args.out}/{g1}')
      for g2, p2 in zip(pheno_2, prefix_2):
        if not os.path.isdir(f'{args.out}/{g1}/{g2}'): os.mkdir(f'{args.out}/{g1}/{g2}')
        if args.corresponding:
            if p2.find(g1.replace('local','')) == -1: continue
            if fnmatch(p2, '*_l_*') or fnmatch(p2, '*_r_*'): continue
        out_rg = f'{args.out}/{g1}/{g2}/{g1}_{p1}.{g2}_{p2}.rg.log'
        
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
        
        if os.path.isfile(out_rg) and (not args.force): continue
        
        submitter.add('bash '+
          f'{scripts_path}/ldsc_master.sh ldsc.py --ref-ld-chr {args.ldsc}/baseline/'+
          f' --w-ld-chr {args.ldsc}/baseline/ '+
          f'--rg {args._in}/{g1}/{p1}.sumstats,'+
          f'{args._in}/{g2}/{p2}.sumstats '+
          f'--out {out_rg[:-4]}')
    
    submitter.submit()
    
if __name__ == '__main__':
    from _utils.slurm import slurm_parser
    parser = slurm_parser(description = 
      'This programme estimates genetic cross-correlation for regional phenotypes')
    parser.add_argument('-p1', nargs = '*', help = 'regional phenotypes',
      default=['deg_local','degi_local','degc_local','clu_local','eff_local','mpl_local'])
    parser.add_argument('-p2', nargs = '*', help = 'correlates phenotype groups',
      default = ['disorders','gradients'])
    parser.add_argument('-i','--in', dest = '_in', help = 'GWA file directory',
      default = '../gcorr/ldsc_sumstats/')
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
      default = '../local_corr/')
    parser.add_argument('-c','--corresponding', dest = 'corresponding', 
      help = 'match local phenotype to corresponding global phenotype',
      default = False, action = 'store_true')
    parser.add_argument('--ldsc', dest = 'ldsc', help = 'LDSC executable directory',
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/ldsc/') # intended to be absolute
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
      default = False, action = 'store_true')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','ldsc']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    
    from _utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheng_%pheno_%maf.sumstats', __file__)
    proj.add_output(args.out+'/%pheng/%pheng/%pheng_%pheno_%maf.%pheng_%pheno_%maf.rg.log', __file__)
    try: main(args)
    except: cmdhistory.errlog()