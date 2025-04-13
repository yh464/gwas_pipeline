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
parser$add_argument('--p1', help = 'Exposure', nargs = '*')
parser$add_argument('--p2', help = 'Outcome') # only one outcome
parser$add_argument('--pval', help = 'P-value threshold', type = 'numeric', default = 5e-8)
parser$add_argument('-o','--out', help = 'Output prefix')
parser$add_argument('-f','--force',dest = 'force', help = 'force overwrite',
                    default = F, action = 'store_true')
args = parser$parse_args(commandArgs(TRUE))
args$out = normalizePath(args$out)
print('Input options')
print(args)

main = function(args){
  if (! dir.exists(dirname(args$out))) dir.create(dirname(args$out))
  
  library(tidyverse)
  library(TwoSampleMR)
  library(httr)
  
  #### pre-process sumstats ####
  # read input extracted SNPs
  gwa = list()
  for (i in c(1:length(args$gwa))) gwa[[i]] = read.delim(args$gwa[i])
  gwa = bind_rows(gwa)
  
  # select SNPs based on outcome
  snps = read.delim(args$clump, sep = '\\s+')$SNP
  gwa = gwa %>% subset(SNP %in% snps)
  
  # harmonise data
  exp = list()
  n_exp = length(args$p1)
  samplesize = list()
  for (i in c(1:length(args$p1))){
    pheno = strsplit(args$p1[i],'/', fixed = T)
    exp[[i]] = gwa %>% subset((Phenotype == pheno[1]) & (Group == pheno[0])) %>%
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
  
  pheno = strsplit(args$p2,'/')
  out = gwa %>% subset((Phenotype == pheno[1]) & (Group == pheno[0])) %>%
    mutate(outcome = args$p2, id.outcome = args$p2) %>% rename(
      effect_allele.outcome = 'A1', other_allele.outcome = 'A2', eaf.outcome = 'AF1',
      beta.outcome = 'BETA', se.outcome = 'SE', pval.outcome = 'P',
      samplesize.outcome = 'N'
    )
  samplesize[[n_exp + 1]]=out$samplesize.outcome%>% max(na.rm = T)
  out = out[,c('SNP','outcome','id.outcome','samplesize.outcome',
               'effect_allele.outcome','other_allele.outcome',
               'eaf.outcome','beta.outcome','se.outcome','pval.outcome')]
  dat = mv_harmonise_data(exp, out)
  
  #### MVMR-Horse ####
  results = paste0(args$out,'_mvmrhorse.txt')
  if (file.exists(results) & !args$force){
    # format data for MVMR_Horse
    betax = dat$exposure_beta; colnames(betax) = sprintf("betaX%i", 1:n_exp)
    betaxse = dat$exposure_se; colnames(betaxse) = sprintf('betaX%ise', 1:n_exp)
    dat_mrhorse = bind_rows(betax, betaxse)
    dat_mrhorse$betaY = dat$outcome_beta; dat_mrhorse$betaYse = dat$outcome_se
    # names of phenotypes are dat$expname$exposure and dat$outname$outcome
    dat$outname = dat$outname %>% separate_wider_delim(outcome, delim='/', names = c('group','phenotype'))
    
    # source MR_horse
    source('https://github.com/aj-grant/mrhorse/raw/refs/heads/main/mr_horse.R')
    res_mrhorse = mvmr_horse(dat_mrhorse)
    # find p-value
    p1tail = res_mrhorse$MR_Coda %>% lapply(as_tibble) %>% bind_rows()
    p1tail = colMeans(p1tail[,-1]>0)
    p1tail[p1tail > 0.5] = 1 - p1tail[p1tail>0.5]
    out_mrhorse = res_mrhorse$MR_Estimate[,-c('Parameter')] %>% add_column(
      Group = dat$outname$group, Phenotype = dat$outname$phenotype, .before = 1
    ) %>% rename(Beta = 'Estimate', SE = 'SD', 
                 CI95_L = '2.5% quantile', CI95_R = '97.5% quantile')
    out_mrhorse$p = p1tail * 2
    out_mrhorse = out_mrhorse %>% mutate(q = p.adjust(p))
    write_delim(out_mrhorse, results, delim='\t', quote = 'needed')
  }
}