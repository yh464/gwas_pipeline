#### Information ####
# A wrapper for Mendelian Randomisation using Causal Analysis Using 
# Summary Effect estimates (CAUSE)
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2024-10-30

#### Reading command line input ####
library(optparse)
library(here) # for portability
optlist = list(
  # input options
  # the filtering for p-values should be moved into the 'batch' file
  make_option('--c1', dest = 'clump1', help = '*CLUMPED* IDP summary stats, one file only'),
  make_option('--g1', dest = 'gwa1', help = 'raw IDP summary stats, one file only'),
  make_option('--c2', dest = 'clump2', help = '*CLUMPED* disorder summary stats, one file only'),
  make_option('--g2', dest = 'gwa2', help = 'raw disorder summary stats, one file only'),
  # does not require the sample size N1 or N2
  
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

#### Input file processing utility ####
merge_gwa_clump = function(gwa, clump, prefix) {
  # gwa file should have 'SNP' column, clump should be the standard PLINK output
  gwa = read.table(gwa, header = T)
  gwa = subset(gwa, gwa$AF1 > 0.01 & gwa$AF1 < 0.99)
  if ('OR' %in% colnames(gwa)) {gwa['BETA'] = log(gwa['OR'])}
  gwa$Phenotype = prefix
  snp = read.table(clump, header = T)['SNP']
  clumped = merge(gwa, snp, by = 'SNP')
  return(list(gwa,clumped))
}

#### main MR execution block ####
main = function(args){
  tic = proc.time()
  #### Required packages ####
  if (! require(tidyverse)){
    install.packages('tidyverse', repos = 'https://cloud.r-project.org')
    library(tidyverse)
  }
  require(ggplot2) # apparently there is a problem with the ggplot2 function 'scale_linewidth_manual'
  library(cause) # throws an error if not installed
  library(data.table) # built-in in R, so this ensures correct R version
  
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
    prefix2,'_cause_cache.rdata'
  )
  
  if (! file.exists(cache_file) | args$force){
    # load input file and merge SNPs
    stat1 = merge_gwa_clump(args$gwa1, args$clump1, prefix1)
    stat2 = merge_gwa_clump(args$gwa2, args$clump2, prefix2)
    stat_fwd = gwas_merge(stat1[[1]], stat2[[1]], snp_name_cols = c('SNP', 'SNP'),
      beta_hat_cols = c('BETA', 'BETA'), se_cols = c('SE', 'SE'), 
      A1_cols = c('A1', 'A1'), A2_cols = c('A2', 'A2'))
    stat_rev = stat_fwd
    colnames(stat_rev) = c('snp','beta_hat_2', 'seb2', 'p2', 'beta_hat_1',
                           'seb1', 'A1', 'A2', 'p1')
    stat_rev = stat_rev[,c(1,5,6,9,2,3,7,8,4)]
    snp_fwd = stat1[[2]]
    snp_fwd = snp_fwd[snp_fwd$SNP %in% stat_fwd$snp,]
    snp_fwd = snp_fwd$SNP
    snp_rev = stat2[[2]]
    snp_rev = snp_rev[snp_rev$SNP %in% stat_fwd$snp,]
    snp_rev = snp_rev$SNP
    
    # estimate CAUSE params
    set.seed(2024)
    varlist = with(stat_fwd, sample(snp, size = 1e+6, replace = F))
    params_fwd = est_cause_params(stat_fwd, varlist)
    # params_rev = est_cause_params(stat_rev, varlist)
    params_rev = params_fwd
    
    # prune SNPs
    stat_fwd = stat_fwd[stat_fwd$snp %in% snp_fwd,]
    stat_rev = stat_rev[stat_rev$snp %in% snp_rev,]
    
    save('stat_fwd', 'stat_rev', 'snp_fwd','snp_rev','params_fwd', 'params_rev','args', file = cache_file)
  } else {load(cache_file)}
  
  toc = proc.time() - tic
  cat('Finished input processing, time = ', toc, 'seconds')
  
  #### Execute CAUSE and write output ####
  # progress check
  out_fwd = paste0(out_prefix,'_mr_forward_cause_results.txt')
  out_rev = paste0(out_prefix,'_mr_reverse_cause_results.txt')
  if (! file.exists(out_fwd) | ! file.exists(out_rev) | args$force) {
    # fit CAUSE
    res_fwd = cause(X = stat_fwd, variants = snp_fwd, param_ests = params_fwd, force = T)
    toc = proc.time() - tic
    cat('Fitted forward CAUSE, time = ', toc, 'seconds')
    
    res_rev = cause(X = stat_rev, variants = snp_rev, param_ests = params_rev, force = T)
    toc = proc.time() - tic
    cat('Fitted reverse CAUSE, time = ', toc, 'seconds')
    
    # parse results
    elpd_fwd = res_fwd$elpd$delta_elpd[c(1,2)]
    summary_fwd = summary(res_fwd)
    tab_fwd = summary_fwd[4]$tab %>% as_tibble()
    tab_fwd$elpd = elpd_fwd
    tab_fwd$p_value = c(NA, summary_fwd[3]$p)
    tab_fwd %>% write.table(file = out_fwd, sep = '\t')
    
    elpd_rev = res_rev$elpd$delta_elpd[c(1,2)]
    summary_rev = summary(res_rev)
    tab_rev = summary_rev[4]$tab %>% as_tibble()
    tab_rev$elpd = elpd_rev
    tab_rev$p_value = c(NA, summary_rev[3]$p)
    tab_rev %>% write.table(file = out_rev, sep = '\t')
  }
}

main(args)
warnings()