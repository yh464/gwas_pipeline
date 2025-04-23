#!/usr/bin/env Rscript
#### Information ####
# Runs multivariable MR for single outcome and multiple exposures
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2025-03-26

#### Command line input ####
library(argparse)
library(here)
parser = ArgumentParser(description = 'This script runs multivariable MR for a single outcome')
parser$add_argument('-i','--in', dest = 'gwa', help = 'GWAS summary stats', nargs = '*')
parser$add_argument('--p1', help = 'Exposure, format <group>/<pheno>, separated by whitespace', nargs = '*')
parser$add_argument('--p2', help = 'Outcome, format <group>/<pheno>, only one outcome')
parser$add_argument('--rg', help = 'Genetic correlation matrix including exposures and outcomes, long format')
parser$add_argument('--pval', help = 'P-value threshold', type = 'numeric', default = 5e-8)
parser$add_argument('-o','--out', help = 'Output prefix')
parser$add_argument('--submodels', dest = 'sub', default = F, action = 'store_true',
  help = 'Test all methods using sub-models identified from MVMR-cML-SuSIE')
parser$add_argument('-f','--force',dest = 'force', help = 'force overwrite',
  default = F, action = 'store_true')
args = parser$parse_args(commandArgs(TRUE))
args$out = normalizePath(args$out)
print('Input options')
print(args)

#### MVMR-Horse model (source from github but modified) ####
mvmr_horse_model = function() {
  for (i in 1:N){
    by[i] ~ dnorm(mu[i], 1 / (sy[i] * sy[i]))
    mu[i] = inprod(bx0[i, 1:K], theta) + alpha[i]
    bx[i,1:K] ~ dmnorm(bx0[i,1:K], Tx[1:K, ((i-1)*K+1):(i*K)])
    
    kappa[i] = (-rho[i]^2 / (1 - K*rho[i]^2))
    bx0[i,1:K] ~ dmnorm(mx + sx0 * rho[i] * alpha[i] / (phi[i] * tau), A - kappa[i] * B)
    r[i] ~ dbeta(10, 10);T(, 1)
    rho[i] = 2*r[i] -1
    alpha[i] ~ dnorm(0, 1 / (tau * tau * phi[i] * phi[i]))
    phi[i] = a[i] / sqrt(b[i])
    a[i] ~ dnorm(0, 1);T(0, )
    b[i] ~ dgamma(0.5, 0.5)
  }
  c ~ dnorm(0, 1);T(0, )
  d ~ dgamma(0.5, 0.5)
  tau = c / sqrt(d)
  mx ~ dmnorm(rep(0, K), R[,])
  for (k in 1:K){
    vx0[k] ~ dnorm(0, 1);T(0, )
    sx0[k] = sqrt(vx0[k])
    theta[k] ~ dunif(-10, 10)
    for (j in 1:K){
      A[j, k] = ifelse(j==k, 1/vx0[j], 0)
      B[j, k] = 1 / (sx0[j] * sx0[k])
    }
  }
}
mvmr_horse_parallel = function(D, no_ini = 3, variable.names = "theta", 
                               n.iter = 10000, n.burnin = 10000){
  if("theta" %in% variable.names){
    variable.names = variable.names
  } else{
    variable.names = c("theta", variable.names)
  }
  
  p = dim(D)[1]
  K = sum(sapply(1:dim(D)[2], function(j){substr(names(D)[j], 1, 5)=="betaX"}))/2
  
  Bx = D[, sprintf("betaX%i", 1:K)]
  Sx = D[, sprintf("betaX%ise", 1:K)]
  Tx = matrix(nrow = K, ncol = p*K)
  for (j in 1:p){
    Tx[, ((j-1)*K+1):(j*K)] = diag(1 / Sx[j, ]^2)
  }
  jags_fit = jags.parallel(data = list(by = D$betaY, bx = Bx, sy = D$betaYse, Tx = Tx, N = p, K = K, R = diag(K)),
                  parameters.to.save = variable.names,
                  n.chains = no_ini,
                  n.iter = n.burnin + n.iter,
                  n.burnin = n.burnin,
                  model.file = mvmr_horse_model)
  mr.coda = as.mcmc(jags_fit)
  mr_estimate = data.frame("Parameter" = sprintf("theta[%i]", 1:K),
    "Estimate" = round(unname(summary(mr.coda)$statistics[sprintf("theta[%i]", 1:K), 1]), 3),
    "SD" = round(unname(summary(mr.coda)$statistics[sprintf("theta[%i]", 1:K), 2]), 3),
    "2.5% quantile" = round(unname(summary(mr.coda)$quantiles[sprintf("theta[%i]", 1:K), 1]), 3),
    "97.5% quantile" = round(unname(summary(mr.coda)$quantiles[sprintf("theta[%i]", 1:K), 5]), 3),
    "Rhat" = round(unname(gelman.diag(mr.coda)$psrf[sprintf("theta[%i]", 1:K), 1]), 3))
  names(mr_estimate) = c("Parameter", "Estimate", "SD", "2.5% quantile", "97.5% quantile", "Rhat")
  return(list("MR_Estimate" = mr_estimate, "MR_Coda" = mr.coda))
}

#### MVMR-Horse main function ####
mvmrhorse = function(dat, subset = NULL, pval = 5e-8, niter = 10000){
  library(R2jags)
  cat('Conducting MVMR-Horse. Time =', proc.time()[3],'\n')
  # takes TwoSampleMR-harmonised data
  # format data for MVMR_Horse
  n_exp = ncol(dat$exposure_beta)
  if (is.null(subset)) subset = 1:n_exp
  # filter p-values
  px = dat$exposure_pval[,subset] %>% matrixStats::rowMins()
  pfilter = px < pval
  if (sum(pfilter) < 5){pthr = px[order(px)[5]]; pfilter = px <= pthr}
  
  betax = dat$exposure_beta[pfilter,subset]
  colnames(betax) = sprintf("betaX%i", 1:ncol(betax))
  betaxse = dat$exposure_se[pfilter,subset]
  colnames(betaxse) = sprintf('betaX%ise', 1:ncol(betaxse))
  dat_mrhorse = bind_cols(betax %>% as.data.frame(row.names = rownames(betax)), 
                          betaxse %>% as.data.frame(row.names = rownames(betaxse)))
  dat_mrhorse$betaY = dat$outcome_beta[pfilter]
  dat_mrhorse$betaYse = dat$outcome_se[pfilter]
  
  # source MR_horse
  # source code has been modified to run parallel JAGS
  set.seed(114514)
  niter = 10000
  res_mrhorse = mvmr_horse_parallel(dat_mrhorse, n.iter = niter)
  # find p-value
  p1tail = res_mrhorse$MR_Coda %>% lapply(as_tibble) %>% bind_rows()
  p1tail = colMeans(p1tail[,-1]>0)
  p1tail[p1tail > 0.5] = 1 - p1tail[p1tail>0.5]
  out_mrhorse = res_mrhorse$MR_Estimate %>% select(-Parameter) %>% add_column(
    Group = dat$expname$group, Phenotype = dat$expname$phenotype, 
    Outcome = dat$outname$outcome, .before = 1
  ) %>% rename(Beta = 'Estimate', SE = 'SD', 
               CI95_L = '2.5% quantile', CI95_R = '97.5% quantile')
  out_mrhorse$p = p1tail * 2
  out_mrhorse = out_mrhorse %>% mutate(q = p.adjust(p))
  if (any(out_mrhorse$Rhat > 1.1)) warning('Model did not converge')
  return(out_mrhorse)
}

#### MVMR-cML-DP sub-function ####
mvmrcmldp = function(dat, rhomat, samplesize, subset = NULL, pval = 5e-8){
  if (is.null(subset)) subset = 1:ncol(dat$exposure_beta)
  rhomat.subset = rhomat[c(subset,n_exp+1), c(subset,n_exp+1)]
  
  # filter summary stats from harmonised data
  px = dat$exposure_pval[,subset] %>% matrixStats::rowMins()
  pfilter = px < pval
  if (sum(pfilter) < 5){pthr = px[order(px)[5]]; pfilter = px <= pthr}
  samplesize.subset = unlist(samplesize)[subset]
  betax = dat$exposure_beta[pfilter, subset]
  betaxse = dat$exposure_se[pfilter, subset]
  betay = dat$outcome_beta[pfilter] %>% as.matrix()
  betayse = dat$outcome_se[pfilter] %>% as.matrix()
  
  # run MVMR-cML
  siginvl = invcov_mvmr(betaxse, betayse, rhomat.subset)
  mvcmlres = MVmr_cML_DP(betax, betay, betaxse, siginvl, min(samplesize.subset),
                         num_pert = 100, K_vec = 0:(length(betay)-3))
  mvcmlres$BIC_se = MVcML_SdTheta(betax, betay, siginvl, mvcmlres$BIC_theta,
                                  zero_ind = setdiff(1:length(betay), mvcmlres$BIC_invalid))
  
  # output table
  mvcmltbl = dat$expname[subset,] %>% select(group, phenotype) %>% add_column(
    total_effect = mvcmlres$BIC_theta[,1], total_se = mvcmlres$BIC_se,
    direct_effect = mvcmlres$BIC_DP_theta[,1], direct_se = mvcmlres$BIC_DP_se[,1]
  ) %>% mutate(total_p = pnorm(-abs(total_effect/total_se))*2,
               direct_p = pnorm(-abs(direct_effect/direct_se))*2) %>%
    mutate(total_q = p.adjust(total_p), direct_q = p.adjust(direct_p))
  return(mvcmltbl)
}

#### MVMR-cML-SuSIE function ####
mvmrsusie = function(dat, rg, samplesize, pval = 5e-8){
  # Takes TwoSampleMR-harmonised data and long-format genetic correlation matrix
  library(MVMRcMLSuSiE)
  library(MVMRcML)
  library(MRcML)
  susie_results = list()
  
  #### Prepare genetic correlation matrix ####
  n_exp = ncol(dat$exposure_beta)
  rhomat = diag(nrow = n_exp + 1)
  names = c(dat$expname$id.exposure, dat$outname$id.outcome)
  colnames(rhomat) = names; rownames(rhomat) = names
  for (i in 1:nrow(rg)) {
    p1 = paste0(rg$group1[i],'/',rg$pheno1[i])
    p2 = paste0(rg$group2[i],'/',rg$pheno2[i])
    if (p1 %in% names & p2 %in% names & p1 != p2){
      rhomat[p1,p2] = rg$rg[i]; rhomat[p2,p1] = rg$rg[i]
    }
  }
  
  #### MVMR-cML-SuSIE step 1 ####
  cat('Conducting step 1 of MVMR-cML-SuSIE, Time =', proc.time()[3],'\n')
  betax = list(); betaxse = list(); betay = list(); betayse = list()
  ## select at least 5 SNPs for each exposure, filtered at p threshold
  for (col in 1:n_exp){
    px = dat$exposure_pval[,col]
    pfilter = px < pval
    if (sum(pfilter) < 5){pthr = px[order(px)[5]]; pfilter = px <= pthr}
    betax[[col]] = dat$exposure_beta[pfilter,col]
    betaxse[[col]] = dat$exposure_se[pfilter,col]
    betay[[col]] = dat$outcome_beta[pfilter]
    betayse[[col]] = dat$outcome_se[pfilter]
  }
  susie.step1 = mvmr.cml.susie.step1(
    sample.sizes = unlist(samplesize), use.openGWAS = F,
    beta.exposure.ls = betax, se.exposure.ls = betaxse,
    beta.outcome.ls = betay, se.outcome.ls = betayse
  )
  
  #### MVMR-cML-SuSIE step 2 ####
  cat('Conducting step 2 of MVMR-cML-SuSIE, Time =', proc.time()[3],'\n')
  ## subset data based on step 1 results
  susie.step1.subset = which(susie.step1 < 0.05/n_exp)
  susie.nsig = length(susie.step1.subset)
  if (susie.nsig == 0) {
    warning('No significant exposure identified by MVMR-cML-SuSIE step 1')
    cat('No significant exposure identified by MVMR-cML-SuSIE step 1\n')
    return(susie_results)
  } else if (susie.nsig == 1) {
    warning('Only one significant exposure identified by MVMR-cML-SuSIE step 1, consider UVMR')
    cat('Only one significant exposure identified by MVMR-cML-SuSIE step 1, consider UVMR\n')
    return(susie_results)
  }
  ## select at least 5 SNPs across all exposures in the subset
  px = dat$exposure_pval[,susie.step1.subset] %>% matrixStats::rowMins()
  pfilter = px < pval
  if (sum(pfilter) < 5){pthr = px[order(px)[5]]; pfilter = px <= pthr}
  samplesize.subset = unlist(samplesize)
  samplesize.subset = samplesize.subset[c(susie.step1.subset, n_exp+1)]
  susie.step2 = mvmr.cml.susie.step2(
    exposure.ids.subset = dat$expname$id.exposure[susie.step1.subset],
    outcome.id = dat$outname$id.outcome,
    sample.sizes.subset = samplesize.subset,
    beta.exposure.mat = dat$exposure_beta[pfilter, susie.step1.subset],
    se.exposure.mat = dat$exposure_se[pfilter, susie.step1.subset],
    pval.exposure.mat = dat$exposure_pval[pfilter, susie.step1.subset],
    beta.outcome.vec = dat$outcome_beta[pfilter],
    se.outcome.vec = dat$outcome_se[pfilter],
    use.openGWAS = F, cutoff = 1 # already filtered for p-value above
  )
  
  #### MVMR-cML-SuSIE step 3 ####
  cat('Conducting step 3 of MVMR-cML-SuSIE, Time =', proc.time()[3],'\n')
  # run preliminarily to determine the number of clusters
  rhomat.subset = rhomat[c(susie.step1.subset, n_exp+1), c(susie.step1.subset, n_exp+1)]
  susie.step3.prelim = mvmr.cml.susie.step3(susie.step2$mvdat, susie.step2$invalid.idx,
                                            susie.step2$theta.vec, rhomat.subset)
  susie.n_clusters = which(susie.step3.prelim$alpha > 1/susie.nsig, arr.ind=T)[,1] %>% 
    unique() %>% length()
  if (susie.n_clusters == 0) {
    warning('No significant signal cluster identified by MVMR-cML-SuSIE step 2')
    cat('No significant signal cluster identified by MVMR-cML-SuSIE step 2\n')
    return(susie_results)
  }
  susie.step3 = mvmr.cml.susie.step3(susie.step2$mvdat, susie.step2$invalid.idx,
    susie.step2$theta.vec, rhomat.subset, susie.n_clusters)
  
  #### compile all combinations from signal clusters ####
  susie.pipmdl = list()
  susie.pip = list()
  for (i in 1:susie.n_clusters) {
    susie.pipmdl[[paste0('exp',i)]] = susie.step1.subset[which(susie.step3$alpha[i,]>1/susie.nsig)]
    susie.pip[[paste0('exp',i)]] = susie.step3$alpha[i,which(susie.step3$alpha[i,]>1/susie.nsig)]
  }
  susie.pip = susie.pip %>% expand.grid() %>% as.matrix() %>% matrixStats::rowProds()
  susie.pipmdl = susie.pipmdl %>% expand.grid() %>% add_column(pip = susie.pip)
  susie_results$models = susie.pipmdl
  if (susie.n_clusters == 1){
    cat('Only one significant signal cluster identified by MVMR-cML-SuSIE, consider UVMR\n')
    susie.pipmdl$exp1 = dat$expname$exposure[susie.pipmdl$exp1]
    return(susie_results)
  }
  
  #### MVMR-cML for each possible model ####
  cat('Testing',nrow(susie.pipmdl),'models from MVMR-cML-SuSIE with MVMR-cML')
  mvcml = list()
  for (i in 1:nrow(susie.pipmdl)){
    # basic information
    cml.subset = susie.pipmdl[i,] %>% select(-pip) %>% as.numeric()
    mvcmltbl = mvmrcmldp(dat, rhomat, samplesize, cml.subset, pval) %>%
      add_column(model_id = i, model_pip = susie.pipmdl$pip[i])
    
    mvcml[[i]] = mvcmltbl
    cat(i, '/', nrow(susie.pipmdl), 'Time =', proc.time()[3],'\n')
  }
  susie_results$res = bind_rows(mvcml)
  return(susie_results)
}

#### main execution block ####
main = function(args){
  if (! dir.exists(dirname(args$out))) dir.create(dirname(args$out))
  
  library(tidyverse)
  library(TwoSampleMR)
  library(matrixStats)
  library(httr)
  
  #### pre-process sumstats ####
  # read input extracted SNPs and rg matrix
  gwa = list()
  for (i in c(1:length(args$gwa))) gwa[[i]] = read.delim(args$gwa[i])
  gwa = bind_rows(gwa) %>% na.omit()
  rg = read.delim(args$rg)
  
  # filter out SNPs within 1 Mb (from cross clumping)
  snp_ind = gwa[,c('SNP','CHR','POS')] %>% distinct() %>% add_column(keep = T) %>%
    arrange(CHR, POS)
  chr_last = 0; pos_last = 0
  for (i in 1:nrow(snp_ind)){
    if (snp_ind$CHR[i] == chr_last & snp_ind$POS[i] < pos_last + 1e6) {
      snp_ind$keep[i] = F
    } else {chr_last = snp_ind$CHR[i]; pos_last = snp_ind$POS[i]}
  }
  snp_ind = snp_ind %>% subset(keep)
  gwa = gwa %>% subset(SNP %in% unique(snp_ind$SNP))
  
  # format exposure data
  exp = list()
  n_exp = length(args$p1)
  samplesize = list()
  for (i in c(1:n_exp)){
    pheno = strsplit(args$p1[i],'/', fixed = T)[[1]]
    exp[[i]] = gwa %>% subset((Phenotype == pheno[2]) & (Group == pheno[1])) %>%
      mutate(exposure = args$p1[i], id.exposure = args$p1[i])
    samplesize[[i]] = exp[[i]]$N %>% max(na.rm = T)
  }
  exp = bind_rows(exp) %>% rename(
    effect_allele.exposure = 'A1', other_allele.exposure = 'A2',
    beta.exposure = 'BETA', se.exposure = 'SE', pval.exposure = 'P',
    eaf.exposure = 'AF1', samplesize.exposure = 'N'
  )
  exp = exp[,c('SNP','exposure','id.exposure','samplesize.exposure',
               'effect_allele.exposure','other_allele.exposure',
               'eaf.exposure','beta.exposure','se.exposure','pval.exposure')]
  n_exp = exp$id.exposure %>% unique() %>% length()
  
  # only select loci above significance level for at least one exposure
  snp_sig = exp %>% group_by(SNP) %>% summarise(P.min = min(pval.exposure)) %>% 
    subset(P.min < args$pval)
  snp_sig = snp_sig$SNP
  exp = exp %>% subset(SNP %in% snp_sig)
  
  # format outcome data
  pheno = strsplit(args$p2,'/')[[1]]
  out = gwa %>% 
    subset((Phenotype == pheno[2]) & (Group == pheno[1]) & (SNP %in% snp_sig)) %>%
    mutate(outcome = args$p2, id.outcome = args$p2) %>% rename(
      effect_allele.outcome = 'A1', other_allele.outcome = 'A2', 
      eaf.outcome = 'AF1', beta.outcome = 'BETA', se.outcome = 'SE', 
      pval.outcome = 'P', samplesize.outcome = 'N'
    )
  samplesize[[n_exp + 1]]=out$samplesize.outcome%>% max(na.rm = T)
  out = out[,c('SNP','outcome','id.outcome','samplesize.outcome',
               'effect_allele.outcome','other_allele.outcome',
               'eaf.outcome','beta.outcome','se.outcome','pval.outcome')]
  
  # harmonise data
  dat = mv_harmonise_data(exp, out)
  # names of phenotypes are dat$expname$exposure and dat$outname$outcome
  dat$expname = dat$expname %>% separate_wider_delim(exposure, delim = '/', 
    names = c('group','phenotype'), cols_remove = F)
  dat$outname = dat$outname %>% separate_wider_delim(outcome, delim='/', 
    names = c('group','phenotype'), cols_remove = F)
  
  #### MVMR-cML-SuSIE ####
  results = paste0(args$out, '_mvmrcmlsusie.txt')
  models = paste0(args$out, '_mvmrcmlsusie_mdl.txt')
  if (!file.exists(results) | args$force){
    out_mrcml = mvmrsusie(dat, rg, samplesize, pval=args$pval)
    if (! is.null(out_mrcml$res)) write_tsv(
      out_mrcml$res, results)
    if (! is.null(out_mrcml$models)) write_tsv(
      out_mrcml$models, models)
  }
  if (file.exists(models)) susie.pipmdl = read_tsv(models) else susie.pipmdl = tibble()
  
  #### MVMR-IVW ####
  results = paste0(args$out,'_mvmrivw.txt')
  if (!file.exists(results) | args$force){
    cat('Performing MVMR-IVW. Time =', proc.time()[3],'\n')
    out_mrivw = mv_multiple(dat)[['result']] %>% 
      select(-id.exposure, -exposure, -id.outcome) %>% rename(p = 'pval') %>%
      mutate(q = p.adjust(p))
    out_mrivw$nsnp = nrow(dat$exposure_beta) # there is a bug in mv_multiple
    print(out_mrivw)
    write_tsv(out_mrivw, results)
  }
  
  #### MVMR-Horse ####
  results = paste0(args$out,'_mvmrhorse.txt')
  if (!file.exists(results) | args$force){
    out_mrhorse = mvmrhorse(dat, pval=args$pval)
    print(out_mrhorse)
    write_tsv(out_mrhorse, results)
  }
  out_mrhorse = read_tsv(results)
  results_50k = paste0(args$out,'_mvmrhorse_50k.txt')
  if (max(out_mrhorse$Rhat) > 1.05 & !file.exists(results_50k)) {
    print('MR-Horse did not converge, retrying with 50000 iterations')
    out_mrhorse = mvmrhorse(dat, pval = args$pval, niter = 50000)
    print(out_mrhorse)
    write_tsv(out_mrhorse, results_50k)
  }
  
  #### MVMR-Horse for all sub-models from MVMR-cML-DP ####
  results_sub = paste0(args$out,'_mvmrhorse_sub.txt')
  if (args$sub & ncol(susie.pipmdl) > 2 & (!file.exists(results_sub) | args$force)) {
    # MVMR-Horse for all sub-models, better than MVMR-cML-DP in simulation
    cat('Conducting MVMR-Horse for all sub-models from MVMR-cML-SuSIE. Time =', proc.time()[3],'\n')
    tmpdir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/temp/mr_cache'
    out_mrhorse_sub = list()
    if (! dir.exists(tmpdir)) dir.create(tmpdir)
    for (i in 1:nrow(susie.pipmdl)){
      subset = susie.pipmdl[i,] %>% select(-pip) %>% as.numeric()
      subresult = paste0(tmpdir,basename(args$out),'_mrhorse_sub',i,'.txt')
      if (file.exists(subresult)) mvcmltbl = read_tsv(subresult) else {
        subtbl = mvmrhorse(dat, subset, args$pval) %>%
          add_column(model_id = i, model_pip = susie.pipmdl$pip[i])
      }
      out_mrhorse_sub[[i]] = subtbl
      cat(i, '/', nrow(susie.pipmdl), 'Time =', proc.time()[3],'\n')
    }
    out_mrhorse_sub %>% bind_rows %>% write_tsv(results_sub)
  }
}

main(args)