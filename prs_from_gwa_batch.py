#!/usr/bin/env python3
'''
This script generates PRScs from external sumstats
'''

def main(args):
    from fnmatch import fnmatch
    import pandas as pd
    
    # array submitter
    from _utils.slurm import array_submitter
    submitter = array_submitter(
        name = 'prs_from_gwa', n_cpu = 2,
        env = 'gentoolspy',
        timeout = 120, lim = 1)
    scripts_path = os.path.realpath(__file__)
    scripts_path = os.path.dirname(scripts_path)
    
    # parse input
    bed_list = open(args.bed).read().splitlines()
    if not os.path.isdir(args.out): os.mkdir(args.out)
    
    tmpdir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/temp/prs_temp'
    if not os.path.isdir(tmpdir): os.mkdir(tmpdir)
    
    for p in args.pheno:
        if not os.path.isdir(f'{tmpdir}/{p}'): os.mkdir(f'{tmpdir}/{p}')
        for x in os.listdir(f'{args._in}/{p}'):
            if fnmatch(x, '*_X.fastGWA') or not fnmatch(x, '*.fastGWA'): continue
            prefix = x.replace('.fastGWA','')
            print(prefix)
            tmpgwa = f'{tmpdir}/{p}/{prefix}.txt'
            tmpn = f'{tmpdir}/{p}/{prefix}_n.txt'
            
            if not os.path.isfile(tmpgwa) or not os.path.isfile(tmpn) or args.force:
                # format GWAS
                df = pd.read_table(f'{args._in}/{p}/{x}')
                if 'OR' in df.columns:
                    tmpdf = df[['SNP','A1','A2','OR','P']]
                elif 'BETA' in df.columns:
                    tmpdf = df[['SNP','A1','A2','BETA','P']]
                elif 'Z' in df.columns:
                    tmpdf = df[['SNP','A1','A2','Z','P']]
                    tmpdf.columns = ['SNP','A1','A2','BETA','P']
                else: print(f'{p}/{x} missing signed statistics, skipping'); continue
                if not 'N' in df.columns and 'N_CAS' in df.columns and 'N_CON' in df.columns:
                    n = df['N_CAS'].max() + df['N_CON'].max()
                else: n = df['N'].max()
                
                # write cache file
                tmpdf.to_csv(tmpgwa, sep = '\t', index = False)
                with open(tmpn, 'w') as n_file: 
                    print(n, file = n_file)
                    n_file.close()
            else:
                n = open(tmpn).read().splitlines()
                n = int(float(n[0]))
            
            outdir = f'{args.out}/{p}/{prefix}'
            if not os.path.isdir(outdir): os.system(f'mkdir -p {outdir}')
      
            for j in range(22):
                out_prefix = f'{outdir}/{prefix}'
                if os.path.isfile(out_prefix+f'_pst_eff_a1_b0.5_phi{args.phi:.0e}_chr{j+1}.txt') and (not args.force):
                    continue
                submitter.add(f'python {args.prscs}/PRScs.py --ref_dir={args.ref} '+
                          f'--bim_prefix={bed_list[j]} --sst_file={tmpgwa} --n_gwas={n} --out_dir={out_prefix} '+
                          f'--chrom={j+1} --phi={args.phi} --seed 114514')
    submitter.submit()
        
if __name__ == '__main__':
    import argparse
    from _utils.slurm import parser_config
    parser = argparse.ArgumentParser(description = 
      'This script generates PRS by continuous shrinkage from external sumstats')
    parser.add_argument('pheno', help = 'Phenotype groups to generate PRS',
      nargs = '*', default = ['disorders','disorders_subtypes'])
    parser.add_argument('-i','--in', dest = '_in', help = 'input directory',
      default = '../gwa/')
    parser.add_argument('--prscs', dest = 'prscs', help = 'directory of PRSCS executable',
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/PRScs/')
    parser.add_argument('--ref', dest = 'ref', help = 'reference panel',
      default = '../params/ldblk_1kg_eur/')
    parser.add_argument('--bed', dest = 'bed', help = 'list of PLINK binaries for target sample',
      default = '../params/bed_files_ukb.txt')
    parser.add_argument('-o','--out', dest ='out', help = 'output directory',
      default = '../prs/prs_effsize/')
    parser.add_argument('--phi', dest = 'phi', 
      help = 'shrinkage parameter, set to 10**-4 to 10**-2',
      type = float, default = 0.01)
    parser.add_argument('--force','-f', dest = 'force', action = 'store_true',
                        default = False, help = 'force overwrite')
    parser = parser_config(parser)
    args = parser.parse_args()
    import os
    for arg in ['_in','out','prscs','bed','ref']:
        exec(f'args.{arg} = os.path.realpath(args.{arg})')
    args.pheno.sort()
    
    from _utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_output(args.out+'/%pheno/.*', __file__)
    try: main(args)
    except: cmdhistory.errlog()