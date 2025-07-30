import pandas as pd
import numpy as np
import scipy.stats as sts

def validate_columns(colnames, silent = False, **kwargs):
  colnames = [c.upper() for c in colnames] # all in uppercase

  # default column names to be interpreted
  colname_dict = dict(
    chrom = ['CHR','CHR_ID','CHROMOSOME','#CHROM'],
    pos = ['BP','POS','POSITION'],
    snp = ['SNP','ID','SNP_ID','SNP_NAME','VARIANT_ID','RS','RSID','RS_ID','MARKERNAME'],
    pval = ['P','PVAL','PVALUE','P_VALUE','P_VAL'],
    n = ['N','NEFF'],
    nca = ['NCASES','N_CASES','NCAS','N_CAS'],
    nco = ['NCONTROLS','N_CONTROLS','NCON','N_CON'],
    af1 = ['AF1','EAF','MAF','EFFECT_ALLELE_FREQ'],
    a1 = ['A1','EFFECT_ALLELE','EFFECTALLELE','ALLELE1','ALT_ALLELE','ALT','EA'],
    a2 = ['A2','OTHER_ALLELE','OTHERALLELE','ALLELE2','REF_ALLELE','REF','NEA'],
    info = ['INFO','IMPINFO','IMPUTATION_INFO','IMP_R2'],
    zscore = ['Z','ZSCORE','Z_SCORE'],
    # this is to distinghish logistic from linear regression
    beta = ['B','BETA','EFFECT_SIZE','EFFECTSIZE','EFFECT','SIGNED_SUMSTAT'],
    odds_ratio = ['OR','LOR','LOGOR','LOG_OR','LOGODDS','LOG_ODDS'], 
    se = ['SE','STDERR','STANDARD_ERROR','STD_ERR','STDERR_BETA','STDERR_OR','STDERR_LOG_OR','STDERR_LOGOR']
  )
  # add user-specified column names
  for key, val in kwargs.items():
    if key in colname_dict.keys():
      if isinstance(val, str): colname_dict[key].append(val.upper())
      else: colname_dict[key] += [c.upper() for c in val]

  description = dict(
    chrom = 'Chromosome', pos = 'Genomic coordinate', snp = 'SNP ID',
    a1 = 'Effect allele', a2 = 'Non-effect allele', af1 = 'Effect allele frequency',
    n = 'Sample size', nca = 'Number of cases', nco = 'Number of controls', 
    info = 'Imputation quality score (r^2)',
    zscore = 'Z-score', se = 'Standard error',
    beta = 'Beta coefficient from linear regression', 
    odds_ratio = 'Odds ratio (or log odds ratio) from logistic regression', pval = 'P-value'
  )

  out_dict = dict(zip(colname_dict.keys(), [None] * len(colname_dict)))

  for group in colname_dict.keys():
    intersect = np.intersect1d(colname_dict[group], colnames)
    if len(intersect) > 1: 
      raise ValueError(f'More than one column found for {description[group]}: {intersect}')
    if len(intersect) == 1:
      out_dict[group] = intersect[0]
      if not silent: print(f'Interpreting {intersect[0]:15s} as {description[group]}')
  print()
  return out_dict

def _convert_z_p_se(df, z = False, p = False, s = False, **kwargs):
  '''Imputes Z-scores, P-values and standard errors from existing columns'''
  df.columns = [c.upper() for c in df.columns]
  cols = validate_columns(df.columns, silent = True, **kwargs)
  def _impute_beta():
    if cols['beta'] == None and cols['odds_ratio'] == 'OR':
      df['BETA'] = np.log(df['OR']); cols['beta'] = 'BETA'
    elif cols['beta'] == None: cols['beta'] = cols['odds_ratio']
  
  if z and cols['zscore'] == None:
    _impute_beta()
    if cols['pval'] != None:
      df['Z'] = sts.norm.isf(df[cols['pval']] / 2) * np.sign(df[cols['beta']])
      cols['zscore'] = 'Z'
    elif cols['se'] != None and cols['beta'] != None:
      df['Z'] = df[cols['beta']] / df[cols['se']]; cols['zscore'] = 'Z'
    else: raise ValueError('Cannot impute z-score')
    print('Imputed z-scores in Z column')
  
  if s and cols['se'] == None:
    _impute_beta()
    if cols['zscore'] != None and cols['beta'] != None:
      df['SE'] = df[cols['beta']] / df[cols['zscore']]
      cols['se'] = 'SE'
    elif cols['pval'] != None and cols['beta'] != None:
      df['SE'] = abs(df[cols['beta']]) / sts.norm.isf(df[cols['pval']] / 2)
      cols['se'] = 'SE'
    else: raise ValueError('Cannot impute standard error')
    print('Imputed standard errors in SE column')

  if p and cols['pval'] == None:
    if cols['zscore'] != None:
      df['P'] = 2 * sts.norm.sf(abs(df[cols['zscore']]))
      cols['pval'] = 'P'
    elif cols['se'] != None and cols['beta'] != None:
      df['P'] = 2 * sts.norm.sf(abs(df[cols['beta']] / df[cols['se']]))
      cols['pval'] = 'P'
    else: raise ValueError('Cannot impute p-value')
    print('Imputed p-values in P column')
  print()
  return df

def _convert_chrom(df, **kwargs):
  '''Converts chromosome column to numeric format'''
  df.columns = [c.upper() for c in df.columns]
  cols = validate_columns(df.columns, silent = True, **kwargs)
  if cols['chrom'] != None:
    if df[cols['chrom']].dtype == 'object':
      df[cols['chrom']] = df[cols['chrom']].str.replace(
        'X', '23').str.replace('Y', '24').str.replace('XY', '25').str.replace('MT', '26')
      df[cols['chrom']] = pd.to_numeric(df[cols['chrom']], errors='coerce')
    elif df[cols['chrom']].dtype == 'int64': pass  # already numeric
    else: raise ValueError(f'Cannot convert {cols["chrom"]} to numeric format')
  return df

def format_gwas(df, beta = 'BETA', odds_ratio = 'OR', 
    *outcols, out = None, accept_z = False, neff_xform = False, all_beta = False, **kwargs):
  '''
  Formats GWAS summary statistics DataFrame
  Arguments:
    df: either a DataFrame or a path to a file
    beta, odds_ratio: column names for beta coefficient and odds ratio
    outcols: columns for output, DO NOT include BETA or OR
    out: output file path
    accept_z: if True, accepts 'z' column as signed statistic
    neff_xform: if True, transforms number of cases and controls to effective sample size, where
      N_eff = 4 * N_cases * N_controls / (N_cases + N_controls)
    all_beta: 
      if True, converts all odds ratios to log(OR)
      if None, does not convert signed statistics
      if False, converts all log(OR) to odds ratios
    kwargs: manually specify column names, accepts:
      chrom, pos, snp, a1, a2, af1, n, nca, nco, beta, odds_ratio, pval, zscore, se
  '''
  # initial preprocessing
  if not isinstance(df, pd.DataFrame):
    df = pd.read_table(df, sep = '\\s+', low_memory = False)
  df.columns = [c.upper() for c in df.columns]
  if 'NEFF_HALF' in df.columns:
    df['NEFF'] = df['NEFF_HALF'] * 2

  if isinstance(outcols[0],list) or isinstance(outcols[0],tuple): outcols = outcols[0]
  outcols = [c.upper() for c in outcols]
  
  # match input and output columns
  print('Specifying output columns')
  out_cols = validate_columns(outcols, **kwargs)
  # fill out SE, P values and Z scores if requested and not present
  df = _convert_z_p_se(df, out_cols['zscore'] != None, out_cols['pval'] != None, out_cols['se'] != None, **kwargs)
  df = _convert_chrom(df, **kwargs)
  print('Parsing input columns')
  orig_cols = validate_columns(df.columns, **kwargs)
  if orig_cols['a1'] != None: df[orig_cols['a1']] = df[orig_cols['a1']].str.upper()
  if orig_cols['a2'] != None: df[orig_cols['a2']] = df[orig_cols['a2']].str.upper()

  # convert number of cases and controls to effective sample size
  if neff_xform and orig_cols['nca'] != None and orig_cols['nco'] != None:
    if orig_cols['n'] != None:
      df[orig_cols['n']] = 4 * df[orig_cols['nca']] * df[orig_cols['nco']] / (df[orig_cols['nca']] + df[orig_cols['nco']])
    else:
      df['N'] = 4 * df[orig_cols['nca']] * df[orig_cols['nco']] / (df[orig_cols['nca']] + df[orig_cols['nco']])
      orig_cols['n'] = 'N'

  # convert signed statistics
  if orig_cols['odds_ratio'] != None:
    or_median = np.median(df[orig_cols['odds_ratio']])
    if 0.9 < or_median < 1.1:
      print(f'Interpreting {orig_cols["odds_ratio"]} as odds ratio on LINEAR SCALE')
      if all_beta:
        df[beta] = np.log(df[orig_cols['odds_ratio']])
        out_cols['beta'] = beta; orig_cols['beta'] = beta # only beta will be outputted
      else: out_cols['odds_ratio'] = odds_ratio # only OR will be outputted
    elif -0.1 < or_median < 0.1:
      print(f'Interpreting {orig_cols["odds_ratio"]} as odds ratio on LOG SCALE')
      if all_beta == False:
        df[odds_ratio] = np.exp(df[orig_cols['odds_ratio']])
        out_cols['odds_ratio'] = odds_ratio; orig_cols['odds_ratio'] = odds_ratio # only OR will be outputted
      else: out_cols['beta'] = beta; orig_cols['beta'] = orig_cols['odds_ratio'] # only beta will be outputted
    else: raise ValueError(f'Median of odds ratio is {or_median:.3f}, cannot determine if it is in log scale')
  # only consider beta if odds ratio is not found
  elif orig_cols['beta'] != None: out_cols['beta'] = beta # only beta will be outputted
  elif accept_z and orig_cols['zscore'] != None:
    out_cols['beta'] = beta; orig_cols['beta'] = orig_cols['zscore'] # z-score will be outputted as though beta
  else: raise ValueError('Cannot determine signed statistic, please specify beta or odds ratio')

  # extract output columns
  out_cols = []; match_cols = []
  for k, v in out_cols.items():
    if v == None: continue
    if orig_cols[k] == None:
      raise ValueError(f'Column {v} not found in input data, please specify manually')
    out_cols.append(v); match_cols.append(orig_cols[k])
  
  out_df = df[match_cols].rename(columns=dict(zip(match_cols, out_cols)))
  if out != None:
    out_df.to_csv(out, sep = '\t', index = False)
    print(f'Output written to {out}')
  
  return out_df

def filter_gwas(df, pval = 1, maf = 0, indel = True, chrom = None, start = None, end = None, 
  ambig = True, info = 0, samplesize = 0, out = None, **kwargs):
  '''
  Filters GWAS summary statistics DataFrame
  Arguments:
    df: either a DataFrame or a path to a file
    pval: maximum p-value to keep
    indel: if False, filters out indels
    chrom: restrict to this chromosome, if None, no restriction
    start, end: restrict to this genomic region, if None, no restriction
    ambig: if False, filters out strand-ambiguous SNPs (MAF > 0.45)
    info: minimum imputation quality score (r^2) to keep
    samplesize: minimum sample size to keep, as a fraction of the maximum sample size in the dataset
    out: output file path
    **kwargs: manually specify column names, accepts:
      chrom, pos, snp, a1, a2, af1, n, nca, nco, beta, odds_ratio, pval, zscore, se
  '''
  if not isinstance(df, pd.DataFrame):
    df = pd.read_table(df, sep = '\\s+', low_memory = False)
  df.columns = [c.upper() for c in df.columns]
  cols = validate_columns(df.columns, silent = True)

  # filter by p-value
  if pval < 1:
    print(f'Filtering by p-value <= {pval:.3e}')
    if cols['pval'] == None and pval < 1:
      raise ValueError('Cannot filter by p-value without P column')
    df = df[df[cols['pval']] <= pval].reset_index(drop = True)
  
  # filter by chromosome and position
  if not isinstance(chrom, list): chrom = [chrom]
  if 'X' in chrom: chrom.remove('X'); chrom.append(23)
  if 'Y' in chrom: chrom.remove('Y'); chrom.append(24)
  if 'XY' in chrom: chrom.remove('XY'); chrom.append(25)
  if 'MT' in chrom: chrom.remove('MT'); chrom.append(26)
  chrom = [int(c) for c in chrom if c is not None and c != '']
  if chrom != None:
    print(f'Filtering for chromosomes {", ".join(chrom)}')
    if cols['chrom'] == None: raise ValueError('Cannot filter by chromosome without CHR column')
    df = _convert_chrom(df, **kwargs)
    df = df[df[cols['chrom']].isin(chrom)].reset_index(drop = True)

    # start and end positions are ignored if chromosome is not specified
    if start != None and end != None:
      print(f'Filtering for SNPs between chr{chrom}:{start} and chr{chrom}:{end}')
      if cols['pos'] == None:
        raise ValueError('Cannot filter by position without BP or POS column')
      df = df[(df[cols['pos']] >= start) & (df[cols['pos']] <= end)].reset_index(drop = True)
  
  # filter out indels
  # do not filter if alleles are not specified
  if not indel and cols['a1'] != None and cols['a2'] != None:
    print('Filtering out indels')
    atcg = ['A', 'T', 'C', 'G']
    df[cols['a1']] = df[cols['a1']].str.upper()
    df[cols['a2']] = df[cols['a2']].str.upper()
    df = df[df[cols['a1']].isin(atcg) & df[cols['a2']].isin(atcg)].reset_index(drop = True)
  
  # filter for MAF
  # do not filter if maf is not specified
  if cols['af1'] != None:
    if maf > 0: 
      print(f'Filtering for MAF = {maf:.3f}')
      df = df[df[cols['af1']] >= maf & df[cols['af1']] <= 1 - maf].reset_index(drop = True)
    if not ambig: 
      print('Filtering out strand-ambiguous SNPs (MAF > 0.45)')
      df = df[df[cols['af1']] < 0.55 & df[cols['af1']] > 0.45].reset_index(drop = True)
    
  # filter by imputation quality score
  # do not filter if info is not specified
  if info > 0 and cols['info'] != None:
    print(f'Filtering for imputation quality score (r^2) >= {info:.3f}')
    df = df[df[cols['info']] >= info].reset_index(drop = True)

  # filter by sample size
  # do not filter if sample size is not specified
  if samplesize > 0 and cols['n'] != None:
    print(f'Filtering for sample size >= {samplesize} * entire cohort size')
    if samplesize > 1: samplesize = 0; Warning('Sample size should be a fraction')
    elif samplesize > 0.9: Warning('Please check the number of SNPs that survive sample size filtering')
    max_samplesize = df[cols['n']].max()
    df = df[df[cols['n']] >= samplesize*max_samplesize].reset_index(drop = True)

  if out != None:
    df.to_csv(out, sep = '\t', index = False)
    print(f'Output written to {out}')
  
  return df