#### Information ####
# A wrapper for Mendelian Randomisation
# Author: Yuankai He (yh464@cam.ac.uk)
# V1.0:   2024-10-15
# V2.0:   2025-01-20 - added MRlap correction
# Notes:  conducts MR by inverse variance, inverse median, Egger, PRESSO and ML
#         if 'clump' parameters are not supplied, the whole GWA file will be used

#### Reading command line input ####
library(optparse)
library(here) # for portability
optlist = list(
  # input options
  # the filtering for p-values should be moved into the 'batch' file
  make_option('--c1', dest = 'clump1', help = '*CLUMPED* trait 1 summary stats, one file only'),
  make_option('--g1', dest = 'gwa1', help = 'raw trait 1 summary stats, one file only'),
  make_option('--n1', help = 'Sample size of trait 1 summary stats'), # default 54030 for functional
  make_option('--c2', dest = 'clump2', help = '*CLUMPED* trait 2 summary stats, one file only'),
  make_option('--g2', dest = 'gwa2', help = 'raw trait 2 summary stats, one file only'),
  make_option('--n2', help = 'Sample size of trait 2 phenotype summary stats'),
  make_option('--nca', help = '# cases of trait 2 phenotype summary stats'),
  make_option('--nco', help = '# controls of trait 2 phenotype summary stats'),
  
  # for MR-lap based correction
  make_option('--pval', help = 'p-value threshold for MR', type = 'double'),
  make_option('--h21', help = 'LDSC h2 estimate for trait 1', type = 'double'),
  make_option('--h2se1', help = 'LDSC h2 std. err. for trait 1', type = 'double'),
  make_option('--h22', help = 'LDSC h2 estimate for trait 2', type = 'double'),
  make_option('--h2se2', help = 'LDSC h2 std. err. for trait 2', type = 'double'),
  make_option('--rglog', help = 'LDSC genetic correlation log file'),
  
  # output directory
  make_option(c('-o','--out'),dest = 'out', help = 'output directory'),
  
  make_option(c('-f','--force'), dest = 'force', help = 'force overwrite', default = F, action = 'store_true')
  # batch file should contain inputs for: phenotype group (e.g. deg_local), prefix (e.g. V1_LH_ROI_0.01),
  # directory to GWA files (../gwa), directory to clumped file (../clump) and p-value threshold (in :.0e format)
  # output prefix should be also determined from these inputs in the 'batch' file
)

args = parse_args(OptionParser(option_list = optlist))
print('Input options')
print(args)

#### Utility to read summary stats ####
merge_gwa_clump = function(gwa, clump, prefix) {
  # gwa file should have 'SNP' column, clump should be the standard PLINK output
  snp = read.table(clump, header = T)['SNP']
  if (is.null(snp)) stop(paste0('Clump file missing ',clump))
  gwa = read.table(gwa, header = T)
  gwa = subset(gwa, gwa$AF1 > 0.01 & gwa$AF1 < 0.99)
  if ('OR' %in% colnames(gwa)) {gwa['BETA'] = log(gwa['OR'])}
  gwa$Phenotype = prefix
  clumped = merge(gwa, snp, by = 'SNP')
  return(list(gwa,clumped))
}

#### Utility to read rg log files ####
parse_rg_log = function(file) {
  tryCatch({
    library(readr)
    data = read_lines(file, skip_empty_rows = T)
    data = strsplit(data[length(data) - 2],'\\s+')
    data = data[[1]]
    results = list()
    results$lambda = as.numeric(data[11])
    results$lambda_se = as.numeric(data[12])
    return(results)
  }, error = function(e) {
    return(list(lambda = NA, lambda_se = NA))
  }
  )
}

#### main MR execution block ####
main = function(args){
#### Required packages ####
if (! require(tidyverse)){
  install.packages('tidyverse', repos = "https://cloud.r-project.org")
  library(tidyverse)
}
library(ggplot2) # apparently there is a problem with the ggplot2 function 'scale_linewidth_manual'
library(TwoSampleMR) # need devtools to install from github so just raise an error
library(MRPRESSO)

#### Input file processing ####
# file names and directory operations
if (!dir.exists(args$out)) dir.create(args$out)
prefix1 = basename(args$gwa1) %>% gsub('.?fastGWA','',.) %>% gsub('.?txt','',.)
prefix2 = basename(args$gwa2) %>% gsub('.?fastGWA','',.) %>% gsub('.?txt','',.)
out_prefix = paste0(args$out,'/',basename(dirname(args$gwa1)),'_',prefix1,'_',prefix2)

# identify common SNPs, use cache to improve loading speed
tmpdir = here('../temp/mr_cache')
if (! dir.exists(tmpdir)) dir.create(tmpdir)
cache_file = paste0(
  tmpdir,'/',
  basename(dirname(args$gwa1)),'_',
  prefix1,'_',
  basename(dirname(args$gwa2)),'_',
  prefix2,'_cache.rdata'
)
if (! file.exists(cache_file) | args$force){
  print(paste0('Reading file: ', args$gwa1))
  stat1 = merge_gwa_clump(args$gwa1, args$clump1, prefix1)
  print(paste0('Reading file: ', args$gwa2))
  stat2 = merge_gwa_clump(args$gwa2, args$clump2, prefix2)
  ref_snp_fwd = merge(stat1[[2]],stat2[[1]], by = 'SNP') # clumped for phenotype 1
  ref_snp_fwd = ref_snp_fwd$SNP
  ref_snp_fwd = data.frame(SNP=ref_snp_fwd)
  print('SNPs used for forward MR')
  print(ref_snp_fwd$SNP)
  stat1_fwd = merge(stat1[[1]], ref_snp_fwd)
  stat2_fwd = merge(stat2[[1]], ref_snp_fwd)

  ref_snp_rev = merge(stat1[[1]],stat2[[2]], by = 'SNP') # clumped for phenotype 2
  ref_snp_rev = ref_snp_rev$SNP
  ref_snp_rev = data.frame(SNP=ref_snp_rev)
  print('SNPs used for reverse MR')
  print(ref_snp_rev$SNP)
  stat1_rev = merge(stat1[[1]], ref_snp_rev)
  stat2_rev = merge(stat2[[1]], ref_snp_rev)
  
# format data
exp1 = format_data(stat1_fwd, type = 'exposure',
                   snp_col = 'SNP',
                   beta_col = 'BETA',
                   se_col = 'SE',
                   effect_allele_col = 'A1', # NB check and double check
                   other_allele_col = 'A2',
                   eaf_col = 'AF1',
                   pval_col = 'P',
                   min_pval = 1e-200,
                   chr_col = 'CHR',
                   pos_col = 'POS')
out1 = format_data(stat1_rev, type = 'outcome',
                   snp_col = 'SNP',
                   beta_col = 'BETA',
                   se_col = 'SE',
                   effect_allele_col = 'A1', # NB check and double check
                   other_allele_col = 'A2',
                   eaf_col = 'AF1',
                   pval_col = 'P',
                   min_pval = 1e-200,
                   chr_col = 'CHR',
                   pos_col = 'POS')

# SUMMARY STATS FOR DISORDERS NEED MANUAL FIDDLING, consider a separate folder for sumstats
exp2 = format_data(stat2_rev, type = 'exposure',
                   snp_col = 'SNP',
                   beta_col = 'BETA',
                   se_col = 'SE',
                   info_col = 'INFO', # THIS IS DIFFERENT FROM THE ABOVE
                   effect_allele_col = 'A1', # NB check and double check
                   other_allele_col = 'A2',
                   eaf_col = 'AF1',
                   pval_col = 'P',
                   min_pval = 1e-200,
                   chr_col = 'CHR',
                   pos_col = 'POS')
out2 = format_data(stat2_fwd, type = 'outcome',
                   snp_col = 'SNP',
                   beta_col = 'BETA',
                   se_col = 'SE',
                   info_col = 'INFO',
                   effect_allele_col = 'A1', # NB check and double check
                   other_allele_col = 'A2',
                   eaf_col = 'AF1',
                   pval_col = 'P',
                   min_pval = 1e-200,
                   chr_col = 'CHR',
                   pos_col = 'POS')

# harmonise data for 2-sample MR
exp1$samplesize.exposure = as.integer(args$n1)
out2$samplesize.outcome = as.integer(args$n2)
mr_fwd_harm = harmonise_data(exp1, out2, action = 2) # action 2 = default, conservative
# pre-calculate variance for directionality test
if (! is.null(args$nca) & ! is.null(args$nco)){
  r2_fwd = get_r_from_lor(lor = mr_fwd_harm$beta.outcome,
                      af = mr_fwd_harm$eaf.outcome,
                      ncase = as.integer(args$nca),
                      ncontrol = as.integer(args$nco),
                      prevalence = 0.0185, # TODO: CHANGE BY DISEASE
                      model = 'logit',
                      correction = T)
} else {
  r2_fwd = get_r_from_pn(mr_fwd_harm$pval.outcome, mr_fwd_harm$samplesize.outcome)
}
r1_fwd = get_r_from_pn(mr_fwd_harm$pval.exposure, mr_fwd_harm$samplesize.exposure)
mr_fwd_harm$r.outcome = r2_fwd
mr_fwd_harm$r.exposure = r1_fwd

# harmonise data for reverse direction
exp2$samplesize.exposure = as.integer(args$n2)
out1$samplesize.outcome = as.integer(args$n1)
mr_rev_harm = harmonise_data(exp2, out1, action = 2)
# pre-calculate variance for directionality test
if (! is.null(args$nca) & ! is.null(args$nco)){
  r2_rev = get_r_from_lor(lor = mr_rev_harm$beta.outcome,
                           af = mr_rev_harm$eaf.outcome,
                           ncase = as.integer(args$nca),
                           ncontrol = as.integer(args$nco),
                           prevalence = 0.0185, # TODO: CHANGE BY DISEASE
                           model = 'logit',
                           correction = T)
} else {
  r2_rev = get_r_from_pn(mr_rev_harm$pval.outcome, mr_rev_harm$samplesize.outcome)
}
r1_rev = get_r_from_pn(mr_rev_harm$pval.exposure, mr_rev_harm$samplesize.exposure)
mr_rev_harm$r.outcome = r1_rev
mr_rev_harm$r.exposure = r2_rev

save('mr_fwd_harm', 'mr_rev_harm','args', file = cache_file)
} else {
  load(cache_file)
}

toc = proc.time()
print(paste0('Finished input data harmonisation, time = ',toc[3]))

#### Define function for MR ####
all_mr_results = function(harm, prefix, ldsc_params) {
  # harm = harmonised data
  # prefix = file name in export, INCLUDING directory
  # ldsc_params is for MR-lap correction
  
  #### perform tests: MR, sensitivity analysis, pleiotropy, leave-one-out, Steiger ####
  # excluding MR-presso test
  res = mr(harm, method_list=c('mr_ivw', 'mr_weighted_median', 'mr_egger_regression'))
  het = mr_heterogeneity(harm, method_list = c('mr_ivw', 'mr_egger_regression'))
  # mr_weighted_median' don't provide heterogeneity test
  pleio = mr_pleiotropy_test(harm)
  single = mr_singlesnp(harm)
  loo = mr_leaveoneout(harm)
  direc = directionality_test(harm)
  
  #### MRlap correction ####
  source('https://github.com/n-mounier/MRlap/raw/refs/heads/master/R/get_correction.R')
  tryCatch({
    mrlap_res = get_correction(
      IVs = harm %>% select(beta.exposure,se.exposure) %>%
        rename(std_beta.exp = beta.exposure, std_SE.exp = se.exposure),
      lambda = ldsc_params$lambda,
      lambda_se = ldsc_params$lambda_se,
      h2_LDSC = ldsc_params$h2_exp,
      h2_LDSC_se = ldsc_params$h2se_exp,
      alpha_obs = res$b[res$method=='Inverse variance weighted'],
      alpha_obs_se = res$se[res$method=='Inverse variance weighted'],
      n_exp = max(harm$samplesize.exposure, na.rm = T) %>% as.numeric(), 
      # as.numeric() required to prevent integer overflow
      n_out = max(harm$samplesize.outcome, na.rm = T) %>% as.numeric(),
      MR_threshold = args$pval,
      verbose = T
    )
    res = res %>% rbind(
      data.frame(
        id.exposure = res$id.exposure[1], id.outcome = res$id.outcome[1],
        outcome = res$outcome[1], exposure = res$exposure[1],
        method = 'MRlap corrected IVW', nsnp = res$nsnp[1],
        b = mrlap_res$alpha_corrected,
        se = mrlap_res$alpha_corrected_se,
        pval = 2*stats::pnorm(-abs(mrlap_res$alpha_corrected/mrlap_res$alpha_corrected_se))
      )
    )
  }, error = function(e) {cat('MRlap correction failed, check for missing data')}
  )
  
  #### outputs for MR ####
  # scatter plot
  scatter = mr_scatter_plot(res,harm)
  ggsave(paste0(prefix,'_scatterplot.pdf'))
  ggsave(paste0(prefix,'_scatterplot.png'))
  remove(scatter)
  
  # forest plot
  forest = mr_forest_plot(single)
  ggsave(paste0(prefix,'_forest.pdf'))
  ggsave(paste0(prefix,'_forest.png'))
  remove(forest)
  
  # leave-one-out plot
  loo_plot = mr_leaveoneout_plot(loo)
  ggsave(paste0(prefix,'_looplot.pdf'))
  ggsave(paste0(prefix,'_looplot.png'))
  remove(loo_plot)
  
  # tabular output
  write.table(res, paste0(prefix,'_results.txt'), sep = '\t')
  
  # tabular outputs for QC tests
  write.table(direc, paste0(prefix,'_dirtest.txt'), sep = '\t')
  write.table(single, paste0(prefix,'_singlesnp.txt'), sep = '\t')
  write.table(pleio, paste0(prefix,'_pleiotropy.txt'), sep = '\t')
  write.table(loo, paste0(prefix,'_lootest.txt'), sep = '\t')
  
  #### perform tests and outputs in MR-presso ####
  # There might be not enough instruments for MR-presso so use a try catch structure
  tryCatch({
    presso_res = run_mr_presso(harm) # use default parameters
    presso_res = presso_res[[1]]
    presso_table = rbind(presso_res$'Main MR results'[,c(2:6)],
                         c('MR-PRESSO Global',NA,NA,
                           presso_res$'MR-PRESSO results'$'Global Test'$RSSobs,
                           presso_res$'MR-PRESSO results'$'Global Test'$Pvalue))
    write.table(presso_table,paste0(prefix,'_presso_results.txt'), sep = '\t')
    
    outlier_log = file(paste0(prefix,'_presso_outlier.txt'), open = 'w')
    writeLines(c('Distortion Coefficient:',
                 as.character(presso_res$'MR-PRESSO results'$'Distortion Test'$'Distortion Coefficient'),
                 'Distortion P-value:',
                 as.character(presso_res$'MR-PRESSO results'$'Distortion Test'$Pvalue),
                 '',
                 'Outliers:',
                 as.character(presso_res$'MR-PRESSO results'$'Distortion Test'$'Outliers Indices')),
               con = outlier_log)
  }, error = function(e) {
    # if an error occurs, write blank files for parsing
    presso_res = matrix(nrow = 3, ncol = 5) %>% as.data.frame()
    colnames(presso_res) = c('MR Analysis','Causal Estimate','Sd','T-stat','P-value')
    write.table(presso_res,paste0(prefix,'_presso_results.txt'), sep = '\t')
    outlier_log = file(paste0(prefix,'_presso_outlier.txt'), open = 'w')
    writeLines(c('Distortion Coefficient:','Distortion P-value:','','','Outliers:',''),
               con = outlier_log)
  })
}

#### Execute MR ####
ldsc_params = parse_rg_log(args$rglog)
ldsc_params_fwd = c(ldsc_params, list(h2_exp = args$h21, h2se_exp = args$h2se1))
ldsc_params_rev = c(ldsc_params, list(h2_exp = args$h22, h2se_exp = args$h2se2))
all_mr_results(mr_fwd_harm,paste0(out_prefix,'_mr_forward'), ldsc_params_fwd)
toc = proc.time()
print(paste0('Finished forward direction MR, time = ',toc[3]))
all_mr_results(mr_rev_harm,paste0(out_prefix,'_mr_reverse'), ldsc_params_rev)
toc = proc.time()
print(paste0('Finished reverse direction MR, time = ',toc[3]))
}

main(args)
warnings() # IMPORTANT as Rscript suppresses warnings