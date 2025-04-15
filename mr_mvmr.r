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
parser$add_argument('-f','--force',dest = 'force', help = 'force overwrite',
                    default = F, action = 'store_true')
args = parser$parse_args(commandArgs(TRUE))
args$out = normalizePath(args$out)
print('Input options')
print(args)

#### MVMR-Horse function ####
mvmrhorse = function(dat){
  # takes TwoSampleMR-harmonised data
  # format data for MVMR_Horse
  n_exp = ncol(dat$exposure_beta)
  betax = dat$exposure_beta; colnames(betax) = sprintf("betaX%i", 1:n_exp)
  betaxse = dat$exposure_se; colnames(betaxse) = sprintf('betaX%ise', 1:n_exp)
  dat_mrhorse = bind_cols(betax %>% as.data.frame(row.names = rownames(betax)), 
                          betaxse %>% as.data.frame(row.names = rownames(betaxse)))
  dat_mrhorse$betaY = dat$outcome_beta; dat_mrhorse$betaYse = dat$outcome_se
  
  # source MR_horse
  source('https://github.com/aj-grant/mrhorse/raw/refs/heads/main/mr_horse.R')
  res_mrhorse = mvmr_horse(dat_mrhorse)
  # find p-value
  p1tail = res_mrhorse$MR_Coda %>% lapply(as_tibble) %>% bind_rows()
  p1tail = colMeans(p1tail[,-1]>0)
  p1tail[p1tail > 0.5] = 1 - p1tail[p1tail>0.5]
  out_mrhorse = res_mrhorse$MR_Estimate %>% select(-Parameter) %>% add_column(
    Group = dat$outname$group, Phenotype = dat$outname$phenotype, .before = 1
  ) %>% rename(Beta = 'Estimate', SE = 'SD', 
               CI95_L = '2.5% quantile', CI95_R = '97.5% quantile')
  out_mrhorse$p = p1tail * 2
  out_mrhorse = out_mrhorse %>% mutate(q = p.adjust(p))
  return(out_mrhorse)
}

#### MVMR-cML-SuSIE function ####
mvmrsusie = function(dat, rg){
  # Takes TwoSampleMR-harmonised data and long-format genetic correlation matrix
  library(MVMRcMLSuSiE)
  library(MVMRcML)
  library(MRcML)
  #### Prepare genetic correlation matrix ####
  n_exp = ncol(dat$exposure_beta)
  rhomat = diag(nrow = n_exp + 1)
  names = c(dat$expname$id.exposure, dat$outname$id.outcome)
  colnames(rhomat) = names; rownames(rhomat = names)
  for (i in 1:nrow(rg)) {
    p1 = paste0(rg$group1[i],'/',rg$pheno1[i])
    p2 = paste0(rg$group2[i],'/',rg$pheno2[i])
    if (p1 %in% names & p2 %in% names & p1 != p2){
      rhomat[p1,p2] = rg$rg[i]; rhomat[p2,p1] = rg$rg[i]
    }
  }
  
  #### MVMR-cML-SuSIE step 1 ####
  cat('Conducting step 1 of MVMR-cML-SuSIE, Time = ', proc.time()[3])
  betax = list(); betaxse = list(); betay = list(); betayse = list()
  ## select at least 5 SNPs for each exposure, filtered at p threshold
  for (col in 1:n_exp){
    px = dat$exposure_pval[,col]
    pfilter = px < args$pval
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
  cat('Conducting step 2 of MVMR-cML-SuSIE, Time = ', proc.time()[3])
  ## subset data based on step 1 results
  susie.step1.subset = which(susie.step1 < 0.05/n_exp)
  susie.nsig = length(susie.step1.subset)
  if (susie.nsig == 0) {
    warning('No significant exposure identified by MVMR-cML-SuSIE step 1')
    return(NULL)
  }
  ## select at least 5 SNPs across all exposures in the subset
  px = dat$exposure_pval[,susie.step1.subset] %>% matrixStats::rowMins()
  pfilter = px < args$pval
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
  cat('Conducting step 3 of MVMR-cML-SuSIE, Time = ', proc.time()[3])
  # run preliminarily to determine the number of clusters
  rhomat.subset = rhomat[c(susie.step1.subset, n_exp+1), c(susie.step1.subset, n_exp+1)]
  susie.step3.prelim = mvmr.cml.susie.step3(susie.step2$mvdat, susie.step2$invalid.idx,
                                            susie.step2$theta.vec, rhomat)
  susie.n_clusters = which(susie.step3.prelim$alpha > 1/susie.nsig, arr.ind=T)[,1] %>% 
    unique() %>% length()
  if (susie.n_clusters == 0) {
    warning('No significant signal cluster identified by MVMR-cML-SuSIE step 2')
    return(NULL)
  }
  susie.step3 = mvmr.cml.susie.step3(susie.step2$mvdat, susie.step2$invalid.idx,
                                     susie.step2$theta.vec, rhomat, susie.n_clusters)
  
  #### compile all combinations from signal clusters ####
  susie.pipmdl = list()
  susie.pip = list()
  for (i in 1:susie.n_clusters) {
    susie.pipmdl[[paste0('exp',i)]] = susie.step1.subset[which(susie.step3$alpha[i,]>1/susie.nsig)]
    susie.pip[[paste0('exp',i)]] = susie.step3$alpha[i,which(susie.step3$alpha[i,]>1/susie.nsig)]
  }
  susie.pip = susie.pip %>% expand.grid() %>% as.matrix() %>% matrixStats::rowProds()
  susie.pipmdl = susie.pipmdl %>% expand.grid() %>% add_column(pip = susie.pip)
  
  #### MVMR-cML for each possible model ####
  cat('Testing',nrow(susie.pipmdl),'models from MVMR-cML-SuSIE with MVMR-cML')
  mvcml = list()
  for (i in 1:nrow(susie.pipmdl)){
    # basic information
    cml.subset = susie.pipmdl[i,] %>% select(-pip) %>% as.numeric()
    rhomat.subset = rhomat[c(cml.subset,n_exp+1), c(cml.subset,n_exp+1)]
    
    # filter summary stats from harmonised data
    px = dat$exposure_pval[,cml.subset] %>% matrixStats::rowMins()
    pfilter = px < args$pval
    if (sum(pfilter) < 5){pthr = px[order(px)[5]]; pfilter = px <= pthr}
    samplesize.subset = unlist(samplesize)[cml.subset]
    betax = dat$exposure_beta[pfilter, cml.subset]
    betaxse = dat$exposure_se[pfilter, cml.subset]
    betay = dat$outcome_beta[pfilter] %>% as.matrix()
    betayse = dat$outcome_se[pfilter] %>% as.matrix()
    
    # run MVMR-cML
    siginvl = invcov_mvmr(betaxse, betayse, rhomat.subset)
    mvcmlres = MVmr_cML_DP(betax, betay, betaxse, siginvl, min(samplesize.subset),
      num_pert = 100, K_vec = 0:(length(betay)-3))
    mvcmlres$BIC_se = MVcML_SdTheta(betax, betay, siginvl, mvcmlres$BIC_theta,
      zero_ind = setdiff(1:length(betay), mvcmlres$BIC_invalid))
    
    # output table
    mvcmltbl = dat$expname[cml.subset,] %>% select(group, phenotype) %>% add_column(
      total_effect = mvcmlres$BIC_theta[,1], total_se = mvcmlres$BIC_se,
      direct_effect = mvcmlres$BIC_DP_theta[,1], direct_se = mvcmlres$BIC_DP_se[,1]
    ) %>% mutate(total_p = pnorm(-abs(total_effect/total_se))*2,
      direct_p = pnorm(-abs(direct_effect/direct_se))*2) %>%
    mutate(total_q = p.adjust(total_p), direct_q = p.adjust(direct_p)) %>%
    add_column(model_id = i, model_pip = susie.pipmdl$pip[i])
    
    mvcml[[i]] = mvcmltbl
    cat(i, '/', nrow(susie.pipmdl), 'Time = ', proc.time()[3])
  }
  mvcml = bind_rows(mvcml)
  return(mvcml)
}

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
  gwa = bind_rows(gwa)
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
  dat$expname = dat$expname %>% separate_wider_delim(exposure, delim = '/', names = c('group','phenotype'))
  dat$outname = dat$outname %>% separate_wider_delim(outcome, delim='/', names = c('group','phenotype'))
  
  #### MVMR-Horse ####
  results = paste0(args$out,'_mvmrhorse.txt')
  if (!file.exists(results) | args$force){
    out_mrhorse = mvmrhorse(dat)
    print(out_mrhorse)
    write_delim(out_mrhorse, results, delim='\t', quote = 'needed')
  }
  
  #### MVMR-cML-SuSIE ####
  results = paste0(args$out, '_mvmrcmlsusie.txt')
  if (!file.exists(results) | args$force){
    out_mrcml = mvmrsusie(dat, rg)
    write_delim(out_mrcml, results, delim = '\t', quote = 'needed')
  }
}

main(args)