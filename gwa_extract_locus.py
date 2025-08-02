#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-07-31

Utility script to extract given locus from GWAS files

Requires following inputs: 
    GWAS summary statistics (single file)
'''

def _extract_locus(df, chr, start, stop, snps = []):
  if len(snps) > 0:
    return df.loc[df.SNP.isin(snps),:]
  else:
    return df.loc[(df.CHR == chr) & (df.POS > start-1) & (df.POS < stop+1),:]
  
def extract_loci(gwa, loci, out):
  '''
  gwa: path to GWAS summary statistics file
  loci: a DataFrame with columns: CHR, START, STOP, or SNPS, or path to a text file with these columns
  out: directory to save extracted loci
  '''
  import pandas as pd
  import os
  from _utils.gwatools import format_gwas
  if not os.path.isdir(out): os.system(f'mkdir -p {out}')
  df = format_gwas(gwa, 'BETA','OR','CHR','SNP','POS','A1','A2','AF1','N','SE','P', silent = True)
  prefix = os.path.basename(gwa).replace('.fastGWA', '')
  if isinstance(loci, str): loci = pd.read_table(loci)
  loci.columns = loci.columns.str.upper()
  for i in range(loci.shape[0]):
    if 'SNPS' in loci.columns and not pd.isna(loci.loc[i, 'SNPS']):
      snps = loci.loc[i, 'SNPS'].split(',')
      locus = _extract_locus(df, None, None, None, snps)
      chrom = locus.CHR.iloc[0]
      start = locus.POS.min(); stop = locus.POS.max()
    else:
      chrom = loci.loc[i, 'CHR']
      start = loci.loc[i, 'START']
      stop = loci.loc[i, 'STOP']
      locus = _extract_locus(df, chrom, start, stop)
    out_prefix = f'{out}/{prefix}_chr{chrom:.0f}_{start:.0f}_{stop:.0f}.txt'
    locus.to_csv(out_prefix, sep = '\t', index = False)
    del locus
    
def main(args = None, **kwargs):
  # only accepts one GWAS file; call array_submitter elsewhere and refer to this script
  import os
  from _utils.gadgets import namespace
  if args == None: args = namespace(**kwargs)
  gwa = f'{args._in}/{args.pheno}.fastGWA'
  out = f'{args.out}/{os.path.dirname(args.pheno)}'
  if os.path.isfile(args.loci):
    extract_loci(gwa, args.loci, out)
  else:
    import pandas as pd
    loci = pd.DataFrame(dict(CHR = args.chr, START = args.start, STOP = args.stop, SNPS = args.snps))
    extract_loci(gwa, loci, out)
  
if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Extract specified loci from GWAS summary statistics')
  parser.add_argument('pheno', help='Phenotype name <group>/<phenotype>')
  parser.add_argument('-i','--in', dest = '_in', help='Input directory containing GWAS files',
    default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/gwa')
  parser.add_argument('-o','--out', required = True, help='Output directory to save extracted loci')
  parser.add_argument('--loci', help='Path to a file with loci definitions or a DataFrame with CHR, START, STOP, SNPS columns')
  parser.add_argument('--chr', type=int, help='Chromosome number')
  parser.add_argument('--start', type=int, help='Start position of the locus')
  parser.add_argument('--stop', type=int, help='End position of the locus')
  parser.add_argument('--snps', nargs='*', default=[], help='List of SNPs to extract')
  args = parser.parse_args()
  
  import os
  args._in = os.path.realpath(args._in); args.out = os.path.realpath(args.out)
  if args.loci != None: args.loci = os.path.realpath(args.loci)

  from _utils import cmdhistory, logger
  logger.splash(args); cmdhistory.log()
  try: main(args)
  except: cmdhistory.errlog()