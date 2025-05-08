#### Information ####
# A wrapper for Mendelian Randomisation
# Author: Yuankai He (yh464@cam.ac.uk)
# V1.0:   2024-10-15
# V2.0:   2025-01-20 - added MRlap correction
# v2.1:   2025-05-07 - merged mr_extract_snp pipeline from MVMR and rg file structure
# Notes:  conducts MR by inverse variance, inverse median, Egger, PRESSO and ML
#         if 'clump' parameters are not supplied, the whole GWA file will be used

#### Reading command line input ####
library(optparse)
library(here) # for portability
optlist = list(
  # input options
  # the filtering for p-values should be moved into the 'batch' file
  make_option('--p1', help = 'Exposure, format <group>/<pheno>'),
  make_option('--p2', help = 'Outcome, format <group>/<pheno>'),
  make_option(c('-i','--in'), dest = 'inst', help = 'Instruments, absolute path separated by colons'),
  make_option('--meta', help = 'Metadata, absolute path separated by colons'),
  make_option('--c1', dest = 'clump1', help = '*CLUMPED* trait 1 summary stats, one file only'),
  make_option('--g1', dest = 'gwa1', help = 'raw trait 1 summary stats, one file only'),
  make_option('--c2', dest = 'clump2', help = '*CLUMPED* trait 2 summary stats, one file only'),
  make_option('--g2', dest = 'gwa2', help = 'raw trait 2 summary stats, one file only'),
  # if metadata is not given, will infer from input instruments file
  
  # for MR-lap based correction
  make_option('--pval', help = 'p-value threshold for MR', type = 'double'),
  make_option('--h21', help = 'LDSC h2 estimate for trait 1', type = 'double'),
  make_option('--h2se1', help = 'LDSC h2 std. err. for trait 1', type = 'double'),
  make_option('--h22', help = 'LDSC h2 estimate for trait 2', type = 'double'),
  make_option('--h2se2', help = 'LDSC h2 std. err. for trait 2', type = 'double'),
  make_option('--gcovint', help = 'LDSC genetic covariance intercept', type = 'double'),
  make_option('--gcintse', help = 'LDSC gcov intercept SE', type = 'double'),
  
  # for MR-APSS based correction
  make_option('--apss', help = 'Apply MR-APSS correction', action = 'store_true', default = F),
  make_option('--ldsc', help = 'LDSC baseline',
    default = here('../toolbox/ldsc/baseline')),
  make_option('--plink', help = 'PLINK executable',
    default = here('../toolbox/plink')),
  make_option('--bfile', help = 'BED file for clumping',
    default = here('../params/fileforclumping_full')),
  
  # output directory
  make_option(c('-o','--out'),dest = 'out', help = 'output directory'),
  
  make_option(c('-f','--force'), dest = 'force', help = 'force overwrite', 
              default = F, action = 'store_true')
)

args = parse_args(OptionParser(option_list = optlist))
print('Input options')
print(args)

#### Parse metadata file ####
parse_metadata = function(file){
  library(tidyverse)
  dat = read_tsv(file)
  colnames(dat) = colnames(dat) %>% tolower()
  if (is.null(dat$group)) dat$group = basename(dirname(file))
  if (is.null(dat$pheno)) dat$pheno = '*' # all phenotypes in the group
  if (is.null(dat$nca)) {dat$nca = NA; dat$nco = NA} else dat$n = dat$nca + dat$nco
  if (is.null(dat$n)) dat$n = NA
  dat = dat %>% select(group, pheno, n, nca, nco)
  return(dat)
}

#### Utility to read summary stats ####
merge_gwa_clump = function(gwa, clump, pheno) {
  library(tidyverse)
  # each gwa file should have 'SNP' column, clump should be the standard PLINK output
  snp = read_tsv(clump)[['SNP']]
  if (is.null(snp)) stop(paste0('Clump file missing ',clump))
  pheno = strsplit(pheno,'/')[[1]]
  if (is.character(gwa)) if (file.exists(gwa)) gwa = read_tsv(gwa)
  gwa = gwa %>% drop_na(BETA, SE) %>% 
    filter(AF1>.01 & AF1<.99 & SNP%in%snp & Group==pheno[1] & Phenotype==pheno[2])
  return(gwa)
}

#### Define function for MR ####
all_mr_results = function(harm, prefix, ldsc_params, apss_params = NULL){
  # harm = harmonised data by TwoSampleMR
  # prefix = output file name, INCLUDING directory
  # ldsc_params is for MR-lap correction
  # apss_params is for MR-APSS correction
  library(ggplot2) # apparently there is a problem with the ggplot2 function 'scale_linewidth_manual'
  library(TwoSampleMR)
  library(MRPRESSO)
  
  #### perform tests: MR, sensitivity analysis, pleiotropy, leave-one-out, Steiger ####
  # excluding MR-presso test
  res = mr(harm, method_list=c('mr_ivw', 'mr_weighted_median', 'mr_egger_regression'))
  het = mr_heterogeneity(harm, method_list = c('mr_ivw', 'mr_egger_regression'))
  # mr_weighted_median' don't provide heterogeneity test
  pleio = mr_pleiotropy_test(harm)
  single = mr_singlesnp(harm)
  loo = mr_leaveoneout(harm)
  direc = directionality_test(harm)
  
  #### if not correct direction, re-test after Steiger filtering ####
  if (! direc$correct_causal_direction) {
    print('Performing Steiger Filtering because directionality test failed')
    harm_filtered = harm %>% steiger_filtering() %>% filter(
      (steiger_dir | steiger_pval > 0.05) & mr_keep
    )
    tryCatch(
      {
        res_filtered = mr_ivw(harm_filtered$beta.exposure,
                              harm_filtered$beta.outcome,
                              harm_filtered$se.exposure,
                              harm_filtered$se.outcome)
      }, error = function(e) {
        res_filtered = list(b = NA, se = NA, pval = NA)
      }
    )
    res = res %>% rbind(data.frame(
      id.exposure = res$id.exposure[1], id.outcome = res$id.outcome[1],
      outcome = res$outcome[1], exposure = res$exposure[1],
      method = 'Steiger-filtered IVW', nsnp = nrow(harm_filtered),
      b = res_filtered$b, se = res_filtered$se, pval = res_filtered$pval
    ))
  }
  
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
  }, error = function(e) {
    cat('MRlap correction failed, check for missing data')
  }
  )
  
  #### MR-APSS correction ####
  if (args$apss) tryCatch({
    library(MRAPSS)
    mrapss_res = MRAPSS(apss_params$dat, exposure = res$exposure[1],
                        outcome = res$outcome[1], C = apss_params$C,
                        Omega = apss_params$Omega, Cor.SelectionBias = T)
    res = res %>% rbind(
      data.frame(
        id.exposure = res$id.exposure[1], id.outcome = res$id.outcome[1],
        outcome = res$outcome[1], exposure = res$exposure[1],
        method = 'MR-APSS', nsnp = nrow(apss_params$dat),
        b = mrapss_res$beta,
        se = mrapss_res$beta.se,
        pval = mrapss_res$pvalue
      )
    )
  }, error = function(e) {
    cat('MR-APSS correction failed, check for missing data\n')
  }
  )
  
  #### Plots for TwoSampleMR ####
  theme_set(theme_classic())
  # scatter plot
  scatter = mr_scatter_plot(res,harm)
  ggsave(paste0(prefix,'_scatterplot.pdf'), width = 4, height = 4)
  remove(scatter)
  
  # forest plot
  forest = mr_forest_plot(single)
  ggsave(paste0(prefix,'_forest.pdf'))
  remove(forest)
  
  # leave-one-out plot
  loo_plot = mr_leaveoneout_plot(loo)
  ggsave(paste0(prefix,'_looplot.pdf'))
  remove(loo_plot)
  
  #### tabular outputs ####
  # tabular output
  write.table(res, paste0(prefix,'_results.txt'), sep = '\t')
  print(res)
  
  # tabular outputs for QC tests
  write.table(het, paste0(prefix, '_heterogeneity.txt'), sep = '\t')
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

#### main MR execution block ####
main = function(args){
  #### Required packages ####
  if (! require(tidyverse)){
    install.packages('tidyverse', repos = "https://cloud.r-project.org")
    library(tidyverse)
  }
  library(ggplot2) # apparently there is a problem with the ggplot2::scale_linewidth_manual
  library(TwoSampleMR) # need devtools to install from github so just raise an error
  library(MRPRESSO)
    
  #### Input file processing ####
  # file names and directory operations
  if (!dir.exists(args$out)) dir.create(args$out)
  p1 = strsplit(args$p1,'/')[[1]]; p2 = strsplit(args$p2,'/')[[1]]
  out_prefix = paste0(args$out,'/',p1[1],'_',p1[2],'_',p2[2])
  
  # identify common SNPs, use cache to improve loading speed
  tmpdir = here('../temp/mr_cache')
  if (! dir.exists(tmpdir)) dir.create(tmpdir)
  cache_file = paste0(tmpdir,'/uvmr_',paste0(args$p1,'.',args$p2,'_cache.rdata') %>%
    gsub('/','_',.))
  
  if (!file.exists(cache_file) | args$force){
    #### Read instruments from mr_extract_snp pipeline ####
    all_instruments = strsplit(args$inst,':')[[1]] %>% lapply(read_tsv) %>% bind_rows()
    gwa1_fwd = merge_gwa_clump(all_instruments, args$clump1, args$p1)
    gwa2_fwd = merge_gwa_clump(all_instruments, args$clump1, args$p2)
    print('SNPs used for forward MR')
    ref_snp_fwd = intersect(gwa1_fwd$SNP, gwa2_fwd$SNP); print(ref_snp_fwd)
    gwa1_fwd = filter(gwa1_fwd, SNP %in% ref_snp_fwd)
    gwa2_fwd = filter(gwa2_fwd, SNP %in% ref_snp_fwd)
    
    gwa1_rev = merge_gwa_clump(all_instruments, args$clump2, args$p1)
    gwa2_rev = merge_gwa_clump(all_instruments, args$clump2, args$p2)
    print('SNPs used for reverse MR')
    ref_snp_rev = intersect(gwa1_rev$SNP, gwa2_rev$SNP); print(ref_snp_rev)
    gwa1_rev = filter(gwa1_rev, SNP %in% ref_snp_rev)
    gwa2_rev = filter(gwa2_rev, SNP %in% ref_snp_rev)
    
    #### format data for TwoSampleMR ####
    exp1 = TwoSampleMR::format_data(gwa1_fwd, type = 'exposure',
      snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
      effect_allele_col = 'A1', # NB check and double check
      other_allele_col = 'A2', eaf_col = 'AF1', pval_col = 'P', 
      min_pval = 1e-200, chr_col = 'CHR', pos_col = 'POS')
    out1 = TwoSampleMR::format_data(gwa1_rev, type = 'outcome',
      snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
      effect_allele_col = 'A1', # NB check and double check
      other_allele_col = 'A2', eaf_col = 'AF1', pval_col = 'P', 
      min_pval = 1e-200, chr_col = 'CHR', pos_col = 'POS')
    
    exp2 = TwoSampleMR::format_data(gwa2_fwd, type = 'exposure',
      snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
      effect_allele_col = 'A1', # NB check and double check
      other_allele_col = 'A2', eaf_col = 'AF1', pval_col = 'P',
      min_pval = 1e-200, chr_col = 'CHR', pos_col = 'POS')
    out2 = TwoSampleMR::format_data(gwa2_rev, type = 'outcome',
      snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
      effect_allele_col = 'A1', # NB check and double check
      other_allele_col = 'A2', eaf_col = 'AF1', pval_col = 'P',
      min_pval = 1e-200, chr_col = 'CHR', pos_col = 'POS')
    
    #### read metadata and get parameters ####
    all_metadata = strsplit(args$meta,':') %>% lapply(parse_metadata) %>% bind_rows()
    metadata = str_split_fixed(c(args$p1, args$p2), '/', 2) %>% as_tibble() %>% 
      setNames(c('group','pheno')) %>% add_column(n=NA, nca=NA, nco=NA)
    for (i in 1:nrow(metadata)){
      group = metadata$group[i]; pheno = metadata$pheno[i]
      match = (all_metadata$group==group) & (all_metadata$pheno %in% c(pheno,'*'))
      if (length(which(match)) > 1) stop('Conflicting metadata found!')
      if (any(match)) metadata[i,3:5] = all_metadata[match,3:5]
    }
    if (is.na(metadata$n[1])) metadata$n[1] = max(gwa1_fwd$N)
    if (is.na(metadata$n[2])) metadata$n[2] = max(gwa2_rev$N)
    
    #### harmonise data for 2-sample MR ####
    getr = function(harm, meta){
      if (!is.na(meta$nca[1]) & !is.na(meta$nco[1])) harm$r.exposure = get_r_from_lor(
        lor = harm$beta.exposure, af = harm$eaf.exposure, ncase = meta$nca[1],
        ncontrol = meta$nco[1], prevalence = meta$nca[1]/(meta$nca[1]+meta$nco[1]),
        model = 'logit', correction = T
      ) else harm$r.outcome = get_r_from_pn(harm$pval.outcome, harm$samplesize.outcome)
      if (!is.na(meta$nca[2]) & !is.na(meta$nco[2])) harm$r.outcome = get_r_from_lor(
        lor = harm$beta.outcome, af = harm$eaf.outcome, ncase = meta$nca[2],
        ncontrol = meta$nco[2], prevalence = meta$nca[2]/(meta$nca[2]+meta$nco[2]),
        model = 'logit', correction = T
      ) else harm$r.outcome = get_r_from_pn(harm$pval.outcome, harm$samplesize.outcome)
      return(harm)
    }
    # harmonise data for forward direction
    exp1$samplesize.exposure = metadata$n[1]
    out2$samplesize.outcome = metadata$n[2]
    mr_fwd_harm = harmonise_data(exp1, out2, action = 2) %>% getr(metadata)
    
    # harmonise data for reverse direction
    exp2$samplesize.exposure = metadata$n[2]
    out1$samplesize.outcome = metadata$n[1]
    mr_rev_harm = harmonise_data(exp2, out1, action = 2) %>% getr(metadata[c(2,1),])
    
    #### pre-calculate parameters for MR-APSS correction ####
    if (args$apss){
      gwa1 = read_tsv(args$gwa1)
      if (!is.null(gwa1$OR)) gwa1 = gwa1 %>% mutate(BETA = log(OR))
      gwa1 = MRAPSS::format_data(gwa1, snp_col = 'SNP', b_col = 'BETA', se_col = 'SE',
        freq_col = 'AF1', A1_col = 'A1', A2_col = 'A2', p_col = 'P', n_col = 'N')
      gwa2 = read_tsv(args$gwa2)
      if (! is.null(gwa2$OR)) gwa2 = gwa2 %>% mutate(BETA = log(OR))
      gwa2 = MRAPSS::format_data(gwa2, snp_col = 'SNP', b_col = 'BETA', se_col = 'SE',
        freq_col = 'AF1', A1_col = 'A1', A2_col = 'A2', p_col = 'P', n_col = 'N')
      apss_params = MRAPSS::est_paras(dat1 = gwa1, dat2 = gwa2, trait1.name = prefix1,
                                      trait2.name = prefix2, ldscore.dir = args$ldsc)
      rm(gwa1); rm(gwa2);
      apss_fwd = list()
      apss_fwd$dat = MRAPSS::clump(apss_params$dat, plink_bin = args$plink, bfile = args$bfile)
      apss_fwd$C = apss_params$C; apss_fwd$Omega = apss_params$Omega
      
      apss_params$dat %>% rename(b.exp = b.out, b.out = b.exp, se.exp = se.out, 
        se.out = se.exp, pval.exp = pval.out, pval.out = pval.exp)
      apss_rev = list()
      apss_rev$dat = MRAPSS::clump(apss_params$dat, plink_bin = args$plink, bfile = args$bfile)
      apss_rev$C = apss_params$C[c(2,1),c(2,1)]
      apss_rev$Omega = apss_params$Omega[c(2,1),c(2,1)]
      # clumping threshold is 1e-5
      rm(apss_params)
      save('mr_fwd_harm', 'mr_rev_harm','apss_fwd', 'apss_rev','args', file = cache_file)
    } else {apss_fwd = NA; apss_rev = NA}
  } else {temp.args = args; load(cache_file); args = temp.args; rm('temp.args')}
  
  toc = proc.time()
  print(paste0('Finished input data harmonisation, time = ',toc[3]))
  
  #### Execute MR ####
  ldsc_params = list(lambda = args$gcovint, lambda_se = args$gcintse)
  ldsc_params_fwd = c(ldsc_params, list(h2_exp = args$h21, h2se_exp = args$h2se1))
  ldsc_params_rev = c(ldsc_params, list(h2_exp = args$h22, h2se_exp = args$h2se2))
  all_mr_results(mr_fwd_harm,paste0(out_prefix,'_mr_forward'), ldsc_params_fwd, apss_fwd)
  toc = proc.time()
  print(paste0('Finished forward direction MR, time = ',toc[3]))
  all_mr_results(mr_rev_harm,paste0(out_prefix,'_mr_reverse'), ldsc_params_rev, apss_rev)
  toc = proc.time()
  print(paste0('Finished reverse direction MR, time = ',toc[3]))
}

main(args)
warnings() # IMPORTANT as Rscript suppresses warnings