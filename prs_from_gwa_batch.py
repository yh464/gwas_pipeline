#!/usr/bin/env python3
'''
This script generates PRScs from external sumstats
'''

def main(args):
  from fnmatch import fnmatch
  import pandas as pd
  
  # parse input
  if not os.path.isdir(args.out): os.mkdir(args.out)
  
  tmpdir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/temp/prs_temp'
  if not os.path.isdir(tmpdir): os.mkdir(tmpdir)
  
  from ._utils.path import find_gwas, find_bed
  pheno = find_gwas(args.pheno, long = True)
  bed_list = find_bed(args.bed)

  # array submitter
  from ._utils.slurm import array_submitter
  submitter = array_submitter(
    name = 'prs_from_gwa_'+pheno[0][0], n_cpu = 2,
    env = 'gentoolspy', timeout = 120)

  for g, p in pheno:
    if not os.path.isdir(f'{tmpdir}/{g}'): os.mkdir(f'{tmpdir}/{g}')
    if os.path.isfile(f'{args._in}/{g}/{p}_noUKBB.fastGWA'): continue # prefer summary stats excluding UKBB
    prefix = p.replace('_noUKBB','')
    print(p)
    tmpgwa = f'{tmpdir}/{g}/{prefix}.txt'
    tmpn = f'{tmpdir}/{g}/{prefix}_n.txt'
    
    if not os.path.isfile(tmpgwa) or not os.path.isfile(tmpn) or args.force:
      # format GWAS
      hdr = open(f'{args._in}/{g}/{p}.fastGWA').readline().replace('\n','').split()
      idx = [hdr.index('SNP'), hdr.index('A1'), hdr.index('A2')]
      if 'OR' in hdr: idx.append(hdr.index('OR'))
      elif 'BETA' in hdr: idx.append(hdr.index('BETA'))
      elif 'Z' in hdr: idx.append(hdr.index('Z'))
      idx.append(hdr.index('P'))
      idx = [x + 1 for x in idx]
      cmd = ['awk', '-v', 'OFS=\'\\t\'', '\'{print'] + [f'${i},' for i in idx[:-1]] + [f'${idx[-1]}'+'}\'', f'{args._in}/{g}/{p}.fastGWA >', tmpgwa]
      os.system(' '.join(cmd))

      df = pd.read_table(f'{args._in}/{g}/{p}.fastGWA', usecols = ['N'] if 'N' in hdr else ['N_CAS','N_CON'])
      if not 'N' in df.columns and 'N_CAS' in df.columns and 'N_CON' in df.columns:
        n = df['N_CAS'].max() + df['N_CON'].max()
      else: n = df['N'].max()
      
      # write cache file
      with open(tmpn, 'w') as n_file: 
        print(n, file = n_file)
        n_file.close()
    else:
      n = open(tmpn).read().splitlines()[0]
    n = int(float(n))

    outdir = f'{args.out}/{g}/{prefix}'
    if not os.path.isdir(outdir): os.system(f'mkdir -p {outdir}')

    for j in range(22):
      out_prefix = f'{outdir}/{prefix}'
      if os.path.isfile(out_prefix+f'_pst_eff_a1_b0.5_phi{args.phi:.0e}_chr{j+1}.txt') and (not args.force):
        continue
      submitter.add(f'python {args.prscs}/PRScs.py --ref_dir={args.ref} '+
        f'--bim_prefix={bed_list[j]} --sst_file={tmpgwa} --n_gwas={int(n)} --out_dir={out_prefix} '+
        f'--chrom={j+1} --phi={args.phi} --seed 114514')
  submitter.submit()
        
if __name__ == '__main__':
    from ._utils.slurm import slurm_parser
    parser = slurm_parser(description = 
      'This script generates PRS by continuous shrinkage from external sumstats')
    parser.add_argument('pheno', help = 'Phenotype groups to generate PRS',
      nargs = '*', default = ['disorders','disorders_subtypes'])
    parser.add_argument('-i','--in', dest = '_in', help = 'input directory',
      default = '../gwa/')
    parser.add_argument('--prscs', dest = 'prscs', help = 'directory of PRSCS executable',
      default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/PRScs/')
    parser.add_argument('--ref', dest = 'ref', help = 'reference panel',
      default = '../params/ref/ldblk/ldblk_1kg_eur/')
    parser.add_argument('--bed', dest = 'bed', help = 'list of PLINK binaries for target sample',
      default = '../params/bed_files_ukb.txt')
    parser.add_argument('-o','--out', dest ='out', help = 'output directory',
      default = '../prs/prs_effsize/')
    parser.add_argument('--phi', dest = 'phi', 
      help = 'shrinkage parameter, set to 10**-4 to 10**-2',
      type = float, default = 0.01)
    parser.add_argument('--force','-f', dest = 'force', action = 'store_true',
                        default = False, help = 'force overwrite')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','prscs','bed','ref']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    args.pheno.sort()
    
    from ._utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_output(args.out+'/%pheno/.*', __file__)
    try: main(args)
    except: cmdhistory.errlog()