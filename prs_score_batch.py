#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2023-07-22
Version 2: 2025-01-23
Version 3: 2026-02-11 *deprecated prs_from_gwa_batch.py, merged into this script*

Batch runs PRS scoring using PRS-cs
'''
import os
import pandas as pd
from _utils import logger
log = logger.logger()

def format_gwas_4prscs(gwa_file, out_prefix):
    if not os.path.isfile(f'{out_prefix}.txt'):
        log.log(f'Formatting GWAS summary statistics for PRScs: {gwa_file}')
        hdr = open(gwa_file).readline().replace('\n','').split()
        idx = [hdr.index('SNP'), hdr.index('A1'), hdr.index('A2')]
        if 'OR' in hdr: idx.append(hdr.index('OR'))
        elif 'BETA' in hdr: idx.append(hdr.index('BETA'))
        elif 'Z' in hdr: idx.append(hdr.index('Z'))
        idx.append(hdr.index('P'))
        idx = [x + 1 for x in idx]
        cmd = ['awk', '-v', 'OFS=\'\\t\'', '\'{print'] + [f'${i},' for i in idx[:-1]] + [f'${idx[-1]}'+'}\'', gwa_file, '>', f'{out_prefix}.txt']
        os.system(' '.join(cmd))
    try:
        n = open(f'{out_prefix}_n.txt').read().splitlines()[0]
    except:
        hdr = open(gwa_file).readline().replace('\n','').split()
        df = pd.read_table(gwa_file, usecols = ['N'] if 'N' in hdr else ['N_CAS','N_CON'])
        if not 'N' in df.columns and 'N_CAS' in df.columns and 'N_CON' in df.columns:
            n = df['N_CAS'].max() + df['N_CON'].max()
        else: n = df['N'].max()
        with open(f'{out_prefix}_n.txt', 'w') as n_file: 
            print(n, file = n_file)
            n_file.close()
    n = int(float(n))
    return n, f'{out_prefix}.txt'

def main(args):
    from _utils.path import find_gwas, find_bed
    from _utils.slurm import array_submitter
    pheno = find_gwas(args.pheno, dirname = args._in, no_ukb = True, long = True)
    submitter = array_submitter(name = 'prs_score'+'_'.join([g for g,_ in pheno]), 
        n_cpu = 2, timeout = 120, env = 'gentoolspy', partition = 'icelake-himem')
    bed_list = find_bed(args.bed)
    ukb_bed_list = find_bed(args.ukb_bed)
    tmpdir = os.path.realpath('../temp/prs_temp'); os.makedirs(tmpdir, exist_ok = True)
    
    for g, p in pheno:
        out_dir = f'{args.out}/{g}/{p}'
        if os.path.realpath(args.bed) == os.path.realpath(args.ukb_bed) or args.use_ukb_effsize:
            # store UKB effect sizes in a special directory
            efs_dir = os.path.realpath(f'../prs/prs_effsize/{g}/{p}')
        else: efs_dir = out_dir
        os.makedirs(out_dir, exist_ok = True)

        # Score by chromosome
        completed = True
        for j in range(22):
            cmds = []
            effsz = f'{efs_dir}/{p}_pst_eff_a1_b0.5_phi{args.phi:.0e}_chr{j+1}.txt'
            if not args.use_ukb_effsize and (not os.path.isfile(effsz) or args.force):
                n, tmpgwa = format_gwas_4prscs(f'{args._in}/{g}/{p}.fastGWA', f'{tmpdir}/{g}/{p}'.replace('_noUKBB',''))
                cmds.append(f'python {args.prscs}/PRScs.py --ref_dir={args.ref} '+
                    f'--bim_prefix={bed_list[j]} --sst_file={tmpgwa} --n_gwas={int(n)} --out_dir={efs_dir}/{p} '+
                    f'--chrom={j+1} --phi={args.phi} --seed 19260817')
            if args.use_ukb_effsize and (not os.path.isfile(effsz) or args.force):
                n, tmpgwa = format_gwas_4prscs(f'{args._in}/{g}/{p}.fastGWA', f'{tmpdir}/{g}/{p}'.replace('_noUKBB',''))
                cmds.append(f'python {args.prscs}/PRScs.py --ref_dir={args.ref} '+
                    f'--bim_prefix={ukb_bed_list[j]} --sst_file={tmpgwa} --n_gwas={int(n)} --out_dir={efs_dir}/{p} '+
                    f'--chrom={j+1} --phi={args.phi} --seed 19260817')
            if os.path.realpath(args.bed) == os.path.realpath(args.ukb_bed) and not os.path.isfile(f'{out_dir}/{p}_pst_eff_a1_b0.5_phi{args.phi:.0e}_chr{j+1}.txt'):
                cmds.append(f'ln -s {efs_dir}/{p}_pst_eff_a1_b0.5_phi{args.phi:.0e}_chr{j+1}.txt {out_dir}/{p}_pst_eff_a1_b0.5_phi{args.phi:.0e}_chr{j+1}.txt')

            out_fname = f'{out_dir}/{p}.chr{j+1}'
            if not os.path.isfile(out_fname+'.sscore') or args.force:
                completed = False
                cmds.append(f'{args.plink} --bfile {bed_list[j]} --chr {j+1} --score {effsz} 2 4 6 center '+
                    f'cols=fid,denom,dosagesum,scoresums --out {out_fname}')
                
            if len(cmds) > 0: 
                completed = False
                submitter.add(*cmds)

        if not completed: continue
        log.log(f'Combining chromosomes for {g}/{p}')
        all_chrs = pd.concat([pd.read_table(f'{out_dir}/{p}.chr{j+1}.sscore', usecols = ['#FID', 'IID', 'SCORE1_SUM']).rename(
            columns = {'#FID': 'FID', 'SCORE1_SUM': f'chr{j+1}'}).sort_values(
            ['FID','IID']).drop_duplicates(['FID','IID']).set_index(['FID','IID']) for j in range(22)], axis = 1).dropna()
        score_total = all_chrs.sum(axis = 1)
        total = pd.DataFrame(dict(score_total = score_total, score_norm = score_total/score_total.std()), index = all_chrs.index)
        total.to_csv(f'{args.out}/{g}/{p}.txt', index = True, sep = '\t')

    submitter.submit()
     
if __name__ == '__main__':
    from _utils.slurm import slurm_parser
    parser = slurm_parser(description = 
      'This script concatenates the PRS-cs output to produce individual level PRS scores')
    parser.add_argument('pheno', help = 'Phenotype groups to generate PRS',
      nargs = '*', default = ['disorders','disorders_subtypes'])
    parser.add_argument('-i','--in', dest = '_in', help = 'input directory',
      default = '../gwa/')
    parser.add_argument('--prscs', dest = 'prscs', help = 'directory of PRSCS executable',
      default = '/home/yh464/rds/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/PRScs/') # intentionally absolute
    parser.add_argument('--plink', dest = 'plink', help = 'Location of plink2 executable',
      default = '/home/yh464/rds/rds-rb643-ukbiobank2/Data_Genetics/plink2') # intentionally absolute
    parser.add_argument('--ref', dest = 'ref', help = 'reference panel',
      default = '/rds/project/rds-Nl99R8pHODQ/ref/ldblk/ldblk_1kg_eur/') # intentionally absolute
    parser.add_argument('--use_ukb_effsize', action = 'store_true', help = 'Use effect sizes estimated using the UKB genotype data')
    parser.add_argument('--bed', dest = 'bed', help = 'PLINK binaries for target sample, list or directory or prefix',
      default = '../params/bed')
    parser.add_argument('--ukb_bed', dest = 'ukb_bed', help = 'PLINK binaries for UKB sample, list or directory or prefix',
      default = '/rds/project/rds-Nl99R8pHODQ/UKB/Imaging_genetics/yh464/bed/')
    parser.add_argument('-o','--out', dest = 'out', help = 'output directory',
      default = '../prs/prs_score/')
    parser.add_argument('--phi', dest = 'phi', help = 'shrinkage parameter used for prscs',
      type = float, default = 0.01)
    parser.add_argument('--force','-f', dest = 'force', action = 'store_true',
                        default = False, help = 'force overwrite')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','plink','bed','ref','prscs']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
    args.pheno.sort()
    
    from _utils import cmdhistory, path, logger
    logger.splash(args)
    cmdhistory.log()
    proj = path.project()
    proj.add_input(args._in+'/%pheno/%pheno*.txt', __file__)
    proj.add_output(args.out+'/%pheno.txt', __file__)
    try: main(args)
    except: cmdhistory.errlog()