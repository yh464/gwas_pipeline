#!/usr/bin/env python3
'''
creates region-to-region genetic correlation matrices for local phenotypes
'''

def main(args):
    import os
    from fnmatch import fnmatch
    
    # array submitter
    from _utils.slurm import array_submitter
    submitter = array_submitter(
        name = 'gcorr_local',
        timeout = 10,mode = 'long')
    
    scripts_path = os.path.realpath(__file__)
    scripts_path = os.path.dirname(scripts_path)
    
    for x in args.pheno:
      ldscdir = f'{args._in}/{x}/'
      gcorrdir = f'{args.out}/{x}/idp_rg/'
      h2dir = f'{args.out}/{x}/h2/'
      global_sumstats = f'{args.glob}/{x}_0.01.sumstats'.replace('local','global')
      
      if not os.path.isdir(gcorrdir): os.system(f'mkdir -p {gcorrdir}')
      if not os.path.isdir(h2dir): os.system(f'mkdir -p {h2dir}')
      
      os.chdir(ldscdir)
      flist = []
      for y in os.listdir():
        if fnmatch(y,'*.sumstats'):
          flist.append(y.replace('.sumstats',''))
      
      for prefix in flist:
        # regional rg w/ global
        grg_fname = f'{gcorrdir}/global.{prefix}.rg'
        if (not os.path.isfile(grg_fname + '.log')) or args.force:
          submitter.add('bash '+
            f'{scripts_path}/ldsc_master.sh ldsc.py --ref-ld-chr {args.ldsc}/baseline/'+
            f' --w-ld-chr {args.ldsc}/baseline/ --rg {ldscdir}/{prefix}.sumstats,{global_sumstats} '+
            f'--out {grg_fname}')
        
        # regional h2
        h2_fname = f'{h2dir}/{prefix}.h2'
        if (not os.path.isfile(h2_fname + '.log')) or args.force:
          if os.path.isfile(f'{ldscdir}/{prefix}.h2.log'):
            os.system(f'cp {ldscdir}/{prefix}.h2.log {h2dir}')
            continue
          submitter.add('bash '+
            f'{scripts_path}/ldsc_master.sh ldsc.py --ref-ld-chr {args.ldsc}/baseline/'+
            f' --w-ld-chr {args.ldsc}/baseline/ --h2 {ldscdir}{prefix}.sumstats '+
            f'--out {h2_fname}') 
      
      # regional correlations
      for i in range(len(flist)):
        for j in range(i):
          skip = False
          if os.path.isfile(f'{gcorrdir}/{flist[i]}.{flist[j]}.rg.log'):
            if os.path.getsize(f'{gcorrdir}/{flist[i]}.{flist[j]}.rg.log') > 2048:
              skip = True
          
          if skip and (not args.force):
            continue
          
          # for each command, print to job file
          submitter.add(
            f'bash {scripts_path}/ldsc_master.sh ldsc.py --ref-ld-chr {args.ldsc}/baseline/'+
            f' --w-ld-chr {args.ldsc}/baseline/ --rg {ldscdir}/{flist[i]}.sumstats,{ldscdir}/{flist[j]}.sumstats '+
            f'--out {gcorrdir}/{flist[i]}.{flist[j]}.rg')
          
    submitter.submit()

if __name__ == '__main__':
    
    import argparse
    from _utils.slurm import parser_config
    
    parser = argparse.ArgumentParser(description = 
      'This programme creates genetic correlation matrices for local phenotypes')
    parser.add_argument('pheno', help = 'Phenotypes', nargs = '*',
      default=['deg_local','degi_local','degc_local','clu_local','eff_local','mpl_local'])
    parser.add_argument('-i','--in', dest = '_in', help = 'summary stats file directory',
      default = '../gcorr/ldsc_sumstats/')
    parser.add_argument('--ldsc', dest = 'ldsc', help = 'LDSC executable directory',
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/ldsc/')
    parser.add_argument('-g','--global', dest = 'glob', help = 'pqath to global summary stats',
      default = '../gcorr/ldsc_sumstats/global/')
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
      default = '../local_corr/')
    parser.add_argument('-f','--force',dest = 'force', help = 'force output',
      default = False, action = 'store_true')
    parser = parser_config(parser)
    args = parser.parse_args()
    
    import os
    for arg in ['_in','out','ldsc','glob']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    
    from _utils import cmdhistory, path
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheng/%reg_%maf.sumstats', __file__)
    proj.add_output(args.out+'/%pheng/rg/%reg.%reg.rg.log', __file__)
    proj.add_output(args.out+'/%pheng/meta/%reg.%pheno.rg.log', __file__)
    proj.add_output(args.out+'/%pheng/h2/%reg.h2.log', __file__)
    try: main(args)
    except: cmdhistory.errlog()