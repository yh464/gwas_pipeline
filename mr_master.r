#### Information ####
# A wrapper for Mendelian Randomisation
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2024-10-15
# Notes:  conducts MR by inverse variance, inverse median, Egger, PRESSO and ML
#         if 'clump' parameters are not supplied, the whole GWA file will be used

#### Reading command line input ####
library(optparse)
library(here) # for portability
optlist = list(
  # input options
  # the filtering for p-values should be moved into the 'batch' file
  make_option('--c1', dest = 'clump1', help = '*CLUMPED* IDP summary stats, one file only'),
  make_option('--g1', dest = 'gwa1', help = 'raw IDP summary stats, one file only'),
  make_option('--n1',dest = 'n1', help = 'Sample size of IDP summary stats'), # default 54030 for functional
  
  make_option('--c2', dest = 'clump2', help = '*CLUMPED* disorder summary stats, one file only'),
  make_option('--g2', dest = 'gwa2', help = 'raw disorder summary stats, one file only'),
  make_option('--n2', dest = 'n2', help = 'Sample size of disorder phenotype summary stats'),
  make_option('--nca', dest = 'nca', help = '# cases of disorder phenotype summary stats'),
  make_option('--nco', dest = 'nco', help = '# controls of disorder phenotype summary stats'),
  
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
  stat1 = merge_gwa_clump(args$gwa1, args$clump1, prefix1)
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

toc = proc.time()-tic
print(paste0('Finished input data harmonisation, time = ',toc[3]))

#### Define function for MR ####
all_mr_results = function(harm, prefix) {
  # harm = harmonised data
  # prefix = file name in export, INCLUDING directory
  
  #### perform tests: MR, sensitivity analysis, pleiotropy, leave-one-out, Steiger ####
  # excluding MR-presso test
  res = mr(harm, method_list=c('mr_ivw', 'mr_weighted_median', 'mr_egger_regression'))
  het = mr_heterogeneity(harm, method_list = c('mr_ivw', 'mr_egger_regression'))
  # mr_weighted_median' don't provide heterogeneity test
  pleio = mr_pleiotropy_test(harm)
  single = mr_singlesnp(harm)
  loo = mr_leaveoneout(harm)
  direc = directionality_test(harm)
  
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
all_mr_results(mr_fwd_harm,paste0(out_prefix,'_mr_forward'))
toc = proc.time()-tic
print(paste0('Finished forward direction MR, time = ',toc[3]))
all_mr_results(mr_rev_harm,paste0(out_prefix,'_mr_reverse'))
toc = proc.time()-tic
print(paste0('Finished reverse direction MR, time = ',toc[3]))
}

main(args)
warnings() # IMPORTANT as Rscript suppresses warnings