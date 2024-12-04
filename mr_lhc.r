#### Information ####
# A wrapper for testing Mendelian Randomisation causality using latent
# heritable confounder (LHC-MR)
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2024-11-04
# Note:   The LHC-MR package is VERY slow and produces a lot of files not prompted (!)
#         need to carefully set WD
#         not really suitable for large-scale analysis, suggest running 'batch' script
#         in 'debug' mode and manually submit interesting ones

#### reading command line input ####
library(optparse)
library(here) # for portability
optlist = list(
  # input options
  make_option('--g1', dest = 'gwa1', help = 'raw IDP summary stats, one file only'),
  make_option('--n1',dest = 'n1', help = 'Sample size of IDP summary stats'), # default 54030 for functional
  make_option('--g2', dest = 'gwa2', help = 'raw disorder summary stats, one file only'),
  make_option('--n2', dest = 'n2', help = 'Sample size of disorder phenotype summary stats'), # defaults to NULL in case N is given as a column
  
  # LDSC scores, needs for Genomic SEM
  make_option('--ldsc', help = 'LD score file (L2), better independent from study cohorts',
    default = here('../params/ldsc_for_gsem/uk10k.l2.ldscore')),
  make_option('--rho',help = 'local LD score file (rho), better independent from study cohorts',
    default = here('../params/ldsc_for_gsem/uk10k.rho.ldscore')),
  make_option('--hm3', help = 'HapMap3 SNP list',
    default = here('../params/ldsc_for_gsem/w_hm3.snplist')),
  make_option('--refld', help = 'reference LD files for ancestry group, directory of 22 gz files',
    default = '/rds/user/yh464/hpc-work/ldsc/baseline/'),
  
  # output directory
  make_option(c('-o','--out'),dest = 'out', help = 'output directory'),
  
  make_option(c('-f','--force'), dest = 'force', help = 'force overwrite', default = F, action = 'store_true')
)

args = parse_args(OptionParser(option_list = optlist))
print('Input options')
print(args)

#### Format GWAS to satisfy LHC requirements ####
format_gwas = function(f, n){
  gwa = read.table(f, header = T)
  gwa = na.omit(gwa)
  if ('OR' %in% colnames(gwa)) {gwa['BETA'] = log(gwa['OR'])
  gwa = subset(gwa, select = -c(OR))}
  colnames(gwa)[colnames(gwa) == 'SE'] = 'SE_BETA'
  if (! 'N' %in% colnames(gwa)) {
    if (! is.null(n)) gwa$N = as.numeric(n) else stop('Sample size for ',f,' is not given')
  }
  return(gwa)
}

#### Main LHC function ####
main = function(args) {
  library(lhcMR)
  library(ieugwasr) # required to remove SNPs not in the EUR reference panel
  if (! require(tidyverse)){
    install.packages('tidyverse', repos = "https://cloud.r-project.org")
    library(tidyverse)
  }
  #### check progress ####
  prefix1 = basename(args$gwa1) %>% gsub('.?fastGWA','',.) %>% gsub('.?txt','',.)
  prefix2 = basename(args$gwa2) %>% gsub('.?fastGWA','',.) %>% gsub('.?txt','',.)
  out_prefix = paste0(args$out,'/',basename(dirname(args$gwa1)),'_',prefix1,'_',prefix2,'_mr_lhc_results.txt')
  if (! file.exists(out_prefix) | args$force){
  
  #### input processing ####
  tic = proc.time()
  tmpdir = here('../temp/mr_cache')
  if (! dir.exists(tmpdir)) dir.create(tmpdir)
  cache = paste0(tmpdir, '/', prefix1,'_',prefix2,'_mr_lhc_cache.rdata')
  setwd(tmpdir) # important because the script automatically writes a fair number of files
  
  if (! file.exists(cache) | args$force) {
    g1 = format_gwas(args$gwa1, args$n1)
    g2 = format_gwas(args$gwa2, args$n2)
    trait.names = c(prefix1, prefix2)
    df = merge_sumstats(list(g1,g2),trait.names,args$ldsc, args$rho)
    # df = df[df$RSID %in% ld_reflookup(df$RSID, opengwas_jwt = 'eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJ5aDQ2NEBjYW0uYWMudWsiLCJpYXQiOjE3MzA4MjQxNDcsImV4cCI6MTczMjAzMzc0N30.ly8hPUsjovjXhWiofO80cpW-mg8ol1SrjVIro8JTegq8P9PW6XiiWI7JNLxiFO3KRQ88Sq5SthUpUG8QV_l5QW1GMvpujlYaI8_QdZkqJa6MO6sgY7BMGW6XYtMjKmMhI_rlM_dPJGbK52A8PYseTb4zsk9BZ22AiIuirZC5xHlBIcqm28Lv1JdrcdeqU2qP9PJIe73urDRUXDSe1ePs_hqJN_7V0g9UQbNiUyG0iAkxZrMr8mtUaL677xAyaBRhpklMpDW73W-Zh22z9jBv1OnCbjLA-RpUJ2WCVUrg66tVX6qcpxsAwj8O-BeGi1l_XU6qwBQIEdMo4gXcwE_U0g'),]
    save(df, args, trait.names, file = cache)
  } else load(cache)
  
  toc = proc.time() - tic
  print(paste0('Finished input processing, time = ', toc[3]))
  
  #### run LHC-MR ####
  if (! file.exists(paste0(tmpdir,'/',prefix1,'_',prefix2,'_mr_lhc_results.rdata')) | args$force){
    sp_list = calculate_SP(df, trait.names, 
                           run_ldsc = F, # FALSE because this will write a whole load of asd2019.sumstats.gz files and not suitable for parallel jobs
                           run_MR = T,
      # there has been a problem with the package where the clump_data function removes all SNPs
      # can be solved by adding an openGWAS JWT to the Renviron file (conda_env/lib/R/etc/Renviron)
                           saveRFiles = F, # we have run MR and LDSC, but these may affect the starting points
                           hm3 = args$hm3, ld = args$refld, nStep = 2,
                           SP_single = 3, SP_pair = 50, SNP_filter = 10)
    toc = proc.time() - tic
    print(paste0('Finished starting point estimation, time = ', toc[3]))
    
    res = lhc_mr(sp_list, trait.names, paral_method = 'lapply', nBlock = 200, nCores = 1) 
    # NB not rslurm because this script is intended to be submitted to SLURM itself
    # Take care when allocating cores (manually set to 4, but may need to change)
    # set nCores to 1 to remove parallelisation, but allocate ~15 GB of memory using 2 CPUs
    # as suggested by the authors
    toc = proc.time() - tic
    print(paste0('Finished MR-LHC, time = ', toc[3]))
    save(res, file = paste0(tmpdir,'/',prefix1,'_',prefix2,'_mr_lhc_results.rdata'))
  } else load(paste0(tmpdir,'/',prefix1,'_',prefix2,'_mr_lhc_results.rdata'))
  write.table(res, file = out_prefix, sep = '\t') 
  # TODO cannot change sep in write.csv, need to manually edit files after current run finishes
  system(paste0('mv ', here(), '/*',prefix1,'* ', tmpdir))
  }
}

main(args)
print(warnings())