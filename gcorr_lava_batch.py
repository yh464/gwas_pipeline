#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-06-12

A flexible framework to conduct regional genetic correlation with LAVA

Requires following inputs: 
  GWAS summary statistics
  LDSC global genetic correlation log
  Clumping output from gwa_clump_batch.py
Outputs:
  local rg, h2 and partial correlation table between one phenotype and all phenotypes of a group
'''

# construct LAVA command for R script
def lava_cmd(out, exp, cov, overlap, clump, task, outfile, 
  tmpdir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/temp/lava/'):
  from hashlib import sha256
  overlap_prefix = sha256((f'{out[0]}:{out[1]})_' + '_'.join([f'{g}:{p}' for g, ps in (exp + cov) for p in ps])
                          ).encode()).hexdigest() + '_overlap.txt'
  overlap_file = f'{tmpdir}/{overlap_prefix}'
  subidx = [f'{out[0]}:{out[1]}'] + [f'{g}:{p}' for g, ps in exp + cov for p in ps]
  overlap_sub = overlap.loc[subidx, subidx]
  overlap_sub.to_csv(overlap_file, sep = '\t', header = True, index = True)
  cmd = ['Rscript', 'gcorr_lava.r', '--p1'] + [f'{g}/{p}' for g, ps in exp for p in ps] + \
        ['--p2', f'{out[0]}/{out[1]}', '--out', outfile, '--overlap', overlap_file] + \
        ['--clump'] + clump + task
  return ' '.join(cmd)

def main(args):
  from _utils.path import find_gwas, find_clump, pair_gwas
  from logparser import crosscorr_parse
  import os
  from statsmodels.stats.moment_helpers import cov2corr
  import pandas as pd
  
  
  # array submitter
  from _utils.slurm import array_submitter
  submitter = array_submitter(name = f'gcorr_lava_{args.p2[0]}_{args.p1[0]}', 
    env = 'gentoolsr', n_cpu = 2, timeout = 240)
  tmpdir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/temp/lava/'
  if not os.path.isdir(tmpdir): os.system(f'mkdir -p {tmpdir}')

  # find files
  exposures = find_gwas(args.p1, dirname = args._in, ext = 'fastGWA', clump = True)
  outcomes = find_gwas(args.p2, dirname = args._in, ext = 'fastGWA', clump = True)
  cov = find_gwas(args.cov, dirname = args._in, ext = 'fastGWA', clump = True)
  meta = []
  for g, _ in exposures + outcomes + cov:
    if os.path.isfile(f'{args._in}/{g}/metadata') and not f'{args._in}/{g}/metadata' in meta: 
      meta.append(f'{args._in}/{g}/metadata')

  # task string general to all commands
  task = ['--ref', args.ref, '--eth', args.eth, '-i', args._in]
  if args.all_loci: task.append('--all-loci')
  if args.all_exp: task.append('--all-exp')
  if args.force: task.append('-f')
  if len(meta) > 0: task += ['--meta'] + meta
  if len(cov) > 0: task += ['--cov'] + [f'{g}:{p}' for g, p in cov]

  if 0 < args.pval < 1 and not args.all_loci:
    exposures_clump = [[find_clump(i, j, dirname = args.clump, pval = args.pval)[0] for j in js] for i, js in exposures]
    outcomes_clump = [[find_clump(i, j, dirname = args.clump, pval = args.pval)[0] for j in js] for i, js in outcomes]
  
  # sample overlap matrix
  gcorr = crosscorr_parse(exposures + outcomes + cov, logdir = args.rg, full = True)
  gcorr['phen1'] = gcorr['group1'] + ':' + gcorr['pheno1']
  gcorr['phen2'] = gcorr['group2'] + ':' + gcorr['pheno2']
  gcorr_rev = gcorr.copy()
  gcorr_rev['phen1'] = gcorr['phen2']; gcorr_rev['phen2'] = gcorr['phen1']
  gcorr = pd.concat([gcorr, gcorr_rev]).drop_duplicates(subset = ['phen1', 'phen2'], keep = 'first')
  gcorr = gcorr.pivot_table(index = 'phen1', columns = 'phen2', values = 'gcov_int')
  overlap = pd.DataFrame(cov2corr(gcorr.values), columns = gcorr.columns, index = gcorr.index)
  overlap.index.name = None

  # all exposures for each outcome
  if args.all_exp:
    for (g2, p2s), c2s in zip(outcomes, outcomes_clump):
      outdir = f'{args.out}/{g2}.'+'_'.join(args.p1)
      if not os.path.isdir(outdir): os.system(f'mkdir -p {outdir}')
      for p2, c2 in zip(p2s, c2s):
        outfile = f'{outdir}/{g2}_{p2}.' + '_'.join(args.p1) + '.lava.txt'
        if os.path.isfile(outfile) and not args.force: continue
        cmd = lava_cmd((g2, p2), exposures, cov, overlap, exposures_clump + [c2], task, outfile)
        submitter.add(cmd)
  else:
    # same file set structure as gcorr_batch.py
    exp_w_clump = [(g1, (p1s, c1s)) for (g1, p1s), c1s in zip(exposures, exposures_clump)]
    out_w_clump = [(g2, (p2s, c2s)) for (g2, p2s), c2s in zip(outcomes, outcomes_clump)]
    pairwise = pair_gwas(exp_w_clump, out_w_clump)

    for g1, (p1s, c1s), g2, (p2s, c2s) in pairwise:
      if g1 > g2: g1, g2, p1s, p2s, c1s, c2s = g2, g1, p2s, p1s, c2s, c1s
      outdir = f'{args.out}/{g1}.{g2}'
      if not os.path.isdir(outdir): os.system(f'mkdir -p {outdir}')
      for p1, c1 in zip(p1s, c1s):
        if g1 == g2:
          p2s = p1s[p1s.index(p1):]; p2s.remove(p1)
          c2s = c1s[p1s.index(p1):]; c2s.remove(c1)
        if len(p2s) == 0: continue
        outfile = f'{outdir}/{g1}_{p1}.{g2}.lava.txt'
        if os.path.isfile(outfile) and not args.force: continue
        # all phenotypes in g2 are included and correlated with g1_p1
        cmd = lava_cmd((g1, p1), [(g2, p2s)], cov, overlap, [c1] + c2s, task, outfile)
        submitter.add(cmd)
  submitter.submit()

if __name__ == '__main__':
  from _utils.slurm import slurm_parser
  parser = slurm_parser(description = 'This script runs LAVA for local genetic correlation')
  parser.add_argument('-p1', nargs = '+', required = True, help = 'Exposures')
  parser.add_argument('-p2', nargs = '*', default = [],
    help = 'Outcomes, leave blank to run pairwise correlations across exposures')
  parser.add_argument('--cov', nargs = '*', default = [], help = 'Covariates, for partial correlation')
  parser.add_argument('-i','--in', dest = '_in', help = 'GWAS summary stats directory',
    default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/gwa')
  parser.add_argument('-o','--out', help = 'Output directory', 
    default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/gcorr/lava/')
  parser.add_argument('--rg', help = 'LDSC genetic correlation log directory',
    default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/gcorr/rglog/')
  parser.add_argument('--clump', dest = 'clump', help = 'Clumping output directory', 
    default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/clump/')
  parser.add_argument('--ref', help = 'Reference LD Blocks directory',
    default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/ref/1000g_by_eth/')
  parser.add_argument('--eth', help = 'Ethnicity', choices = ['eas', 'afr', 'eur', 'sas', 'amr']
    , default = 'eur')
  parser.add_argument('--all-loci', action = 'store_true', help = 'Analyse all loci')
  parser.add_argument('--all-exp', action = 'store_true', help = 'Multiple regression with all exposures')
  parser.add_argument('--pval', type = float, default = 5e-8, help = 'Clumping p-value threshold')
  parser.add_argument('-f', '--force', action = 'store_true', help = 'Force overwrite')
  args = parser.parse_args()

  import os
  for arg in ['_in','out','rg','clump', 'ref']:
    setattr(args, arg, os.path.realpath(getattr(args, arg)))
  
  from _utils import cmdhistory, path, logger
  logger.splash(args)
  cmdhistory.log()
  proj = path.project()
  proj.add_input(f'{args._in}', __file__)
  proj.add_input(f'{args.clump}',__file__)
  proj.add_input(args.rg,__file__)
  proj.add_output(f'{args.out}',__file__)
  
  try: main(args)
  except: cmdhistory.errlog()