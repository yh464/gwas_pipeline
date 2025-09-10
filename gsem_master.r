#!/usr/bin/env Rscript
#### Information ####
# A flexible framework to run genomic SEM
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2025-04-26

#### Command line input ####
library(argparse)
library(here)
parser = ArgumentParser(description = 'This script runs genomic SEM')
# path specs
parser$add_argument('-i','--in', dest = 'input', help = 'input MUNGED summary stats directory',
  default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/gcorr/ldsc_sumstats')
parser$add_argument('--full', dest = 'full', 
  help = 'input FULL summary stats directory, needed for common factor GWAS/GWAS by subtraction',
  default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/gwa')
parser$add_argument('--p1', nargs = '+', 
  help = 'Exposure, format <group>/<pheno>, separated by whitespace')
parser$add_argument('--p2', nargs = '*',
  help = 'Outcome, format <group>/<pheno>, leave blank for common factor analysis')
parser$add_argument('--cov', nargs = '*',
  help = 'Covariates, format <group>/<pheno>')
parser$add_argument('--med', nargs = '*', 
  help = 'Mediators, format <group>/<pheno>')
parser$add_argument('--meta', nargs = '*', help = 'Metadata files')
parser$add_argument('--ld', help = 'LD reference',
  default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/ldsc/baseline')
parser$add_argument('--ref', help = 'Reference file for SNP variance calculation', default = 
  '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/ldsc_for_gsem/ref.1000G.txt')
parser$add_argument('-o','--out', help = 'Output prefix; NB manual models will be output where model file is')

# analyses
parser$add_argument('--common', default = F, action = 'store_true', help = 'Common factor model')
parser$add_argument('--efa', default = F, action = 'store_true', help = 'Exploratory factor analysis')
parser$add_argument('--efa_n', default = -1, type = 'integer', help = 'number of factors')
parser$add_argument('--efa_thr', default = 0.3, help = 'loading threshold to keep an item')
parser$add_argument('--mdl', default = F, action = 'store_true', 
  help = 'Causal model and subtraction model; enter --mdl --gwas for GWAS by subtraction')
parser$add_argument('--manual', help = 'Enter text file containing manually specified model')
parser$add_argument('--gwas', default = 0, type = 'integer', help = 'Output GWAS, by chromosome')
# need to do GWAS by single chromosomes because HPC will time out
parser$add_argument('-f','--force',dest = 'force', help = 'force overwrite',
  default = F, action = 'store_true')
args = parser$parse_args(commandArgs(TRUE))
args$p1 = sort(args$p1); args$p2 = sort(args$p2)
args$cov = sort(args$cov); args$med = sort(args$med)
if (!is.null(args$manual)) args$manual = normalizePath(args$manual)
#### Sanity checks to command line arguments ####
args$input = normalizePath(args$input); args$out = normalizePath(args$out)
if (length(args$meta) < length(c(args$p1,args$cov, args$p2, args$med))) warning(
  'Missing metadata, assuming traits to be continuous')
print('Input options')
print(args)
if (length(args$p2) > 1) stop('Only one or zero outcome phenotypes accepted')
if (!args$mdl & is.null(args$manual) & length(args$p2) == 1) warning(
  'Unneeded p2 parameter supplied')
if (!args$mdl & is.null(args$manual) & length(args$med) > 0) warning(
  'Unneeded med parameter supplied')
if (!args$common & !args$efa & !args$mdl & is.null(args$manual)) stop(
  'Please tell me to do something')

#### paLDSC function to extract number of factors, modified for automation ####
paLDSC_mod = function(S, V) {
  library(MASS)
  library(matrixStats)
  library(gdata)
  library(psych)
  library(Matrix)
  # default parameters
  r = 500; p = 0.95; fm = 'minres'; nfactors = 1;
  
  #### Get Dimensions of Matrices ####
  k=dim(S)[1] #k phenotypes
  kstar=k*(k+1)/2 #kstar unique variances/covariances
  Svec=lowerTriangle(S,diag=T) #vectorize S
  SNULL=(0*S)
  diag(SNULL)=1  
  SNULLvec=lowerTriangle(SNULL,diag=T) #vectorize S null
  Observed_PCA_values=eigen(S)$values
  obs = data.frame(Observed_PCA_values)
  obs$type = c('Observed Data')
  obs$num = c(row.names(obs))
  obs$num = as.numeric(obs$num)
  colnames(obs) = c('eigenvalue', 'type', 'num')
  
  #### Factor parallel analysis ####
  Ssmooth<-as.matrix((nearPD(S, corr = T))$mat) #Smooth S matrix
  Observed_FA_values=fa(Ssmooth, fm = fm, nfactors = nfactors,SMC = FALSE, 
                        warnings = FALSE, rotate = "none")$values
  obsFA = data.frame(Observed_FA_values)
  obsFA$type = c('Observed Data')
  obsFA$num = c(row.names(obsFA))
  obsFA$num = as.numeric(obsFA$num)
  colnames(obsFA) = c('eigenvalue', 'type', 'num')
  EIGfa=as.data.frame(matrix(NA,nrow=k,ncol=r))
  for (i in 1:r) {
    Sample_null=mvrnorm(n=1,mu=SNULLvec,Sigma=V)
    Sample_null_M=matrix(0,ncol=k,nrow=k)
    lowerTriangle(Sample_null_M,diag=T)=Sample_null
    upperTriangle(Sample_null_M,diag=F)=upperTriangle(t(Sample_null_M))
    Ssmooth<-as.matrix((nearPD(Sample_null_M, corr = T))$mat)
    EIGfa[,i]=fa(Ssmooth, fm = fm, nfactors = nfactors,SMC = FALSE, 
                 warnings = FALSE, rotate = "none")$values 
  } 
  
  Parallel_fa_values=rowQuantiles(as.matrix(EIGfa),probs=p)
  simPAfa = data.frame(Parallel_fa_values)
  simPAfa$type = paste("Simulated data (",(p*100),"th %ile)",sep = "")
  simPAfa$num = c(row.names(obs))
  simPAfa$num = as.numeric(simPAfa$num)
  colnames(simPAfa) = c('eigenvalue', 'type', 'num')
  eigendatPAfa = rbind(obsFA,simPAfa)
  nfactPAfa <- min(which((eigendatPAfa[1:k,1] < eigendatPAfa[(k+1):(k*2),1]) == TRUE))-1
  if(nfactPAfa==Inf){
    nfactPAfa <- 1
  }
  obsPAfa <- data.frame(obsFA[1]-simPAfa[1])
  obsPAfa$type = paste("Observed minus simulated data (",(p*100),"th %ile)",sep = "")
  obsPAfa$num = c(row.names(obsFA))
  obsPAfa$num = as.numeric(obsPAfa$num)
  colnames(obsPAfa) = c('eigenvalue', 'type', 'num')
  nfactobsPAfa <- which(obsPAfa < 0)[1]-1
  return(nfactPAfa)
}

#### Parse metadata file ####
parse_metadata = function(file){
  library(tidyverse)
  dat = read_tsv(file)
  colnames(dat) = colnames(dat) %>% tolower()
  if (is.null(dat$group)) dat$group = basename(dirname(file))
  if (is.null(dat$pheno)) dat$pheno = '*' # all phenotypes in the group
  if (is.null(dat$sample_prev)) dat$sample_prev = NA
  if (is.null(dat$pop_prev)) dat$pop_prev = NA
  if (is.null(dat$selogit)) dat$selogit = F else dat$selogit = replace_na(dat$selogit,F)
  if (is.null(dat$ols)) dat$ols = T else dat$ols = replace_na(dat$ols, T)
  if (is.null(dat$linprob)) dat$linprob = F else dat$linprob = replace_na(dat$linprob,F)
  if (is.null(dat$n)) dat$n = NA
  dat = dat %>% select(group, pheno, n, sample_prev, pop_prev, selogit, ols, linprob)
  return(dat)
}

#### Function to write genomic SEM output ####
write_model = function(dwls, file){
  if (is.null(dwls$results)) return(F)
  out = dwls$results %>% rename(
    beta_unstd = 'Unstand_Est', se_unstd = 'Unstand_SE', beta = 'STD_Genotype', 
    se = 'STD_Genotype_SE',  beta_std_full = 'STD_All', p = 'p_value') %>% 
    mutate(p = as.numeric(p) %>% replace_na(0)) %>% mutate(q = p.adjust(p,'BH')) %>% 
    add_column(chi2 = dwls$modelfit$chisq, chi2_p = dwls$modelfit$p_chisq, 
      df = dwls$modelfit$df, CFI = dwls$modelfit$CFI, SRMR = dwls$modelfit$SRMR)
  write_tsv(out, file = file, quote = 'needed')
  return(T)
}

#### main execution block ####
main = function(args){
  library(httr)
  library(GenomicSEM)
  library(openssl)
  library(tidyverse)
  tmpdir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/temp/gsem'
  if (!dir.exists(tmpdir)) dir.create(tmpdir)
  
  #### Read metadata ####
  all_metadata = lapply(args$meta, parse_metadata) %>% bind_rows()
  metadata = str_split_fixed(c(args$p1, args$cov, args$p2, args$med), '/', 2) %>% 
    as_tibble() %>% setNames(c('group','pheno')) %>% 
    add_column(n=NA, sample_prev=NA, pop_prev=NA, selogit=NA, ols=NA, linprob=NA)
  for (i in 1:nrow(metadata)){
    group = metadata$group[i]; pheno = metadata$pheno[i]
    match = (all_metadata$group==group) & (all_metadata$pheno %in% c(pheno,'*'))
    if (length(which(match)) > 1) stop('Conflicting metadata found!')
    if (any(match)) metadata[i,3:8] = all_metadata[match,3:8]
  }
  # defaults to quantitative traits
  metadata$selogit = replace_na(metadata$selogit,F)
  metadata$ols = replace_na(metadata$ols, T)
  metadata$linprob = replace_na(metadata$linprob,F)
  
  #### Prepare modal parameters ####
  # correct trait names for lavaan syntax
  trait.names = gsub('/','_',c(args$p1, args$cov, args$p2))
  trait.names_med = gsub('/','_',c(args$p1, args$cov, args$p2, args$med))
  p1 = gsub('/','_',args$p1); p2 = gsub('/','_',args$p2); 
  cov = gsub('/','_', args$cov); med = gsub('/','_',args$med)
  n = length(c(p1,cov,p2))
  ldsc_cache = paste0(tmpdir,'/',sha256(paste(trait.names_med,collapse='.')),'.ldsc.rdata')
  if (file.exists(ldsc_cache)) load(ldsc_cache) else {
    ldscoutput = ldsc(
      traits = paste0(args$input, '/', c(args$p1, args$cov, args$p2,args$med),'.sumstats'),
      sample.prev = metadata$sample_prev[1:n], population.prev = metadata$pop_prev[1:n],
      ld = args$ld, wld = args$ld, trait.names = trait.names_med,
      ldsc.log = paste0(tmpdir,'/',paste(rand_bytes(4), collapse = '')))
    save(ldscoutput, file = ldsc_cache)
  }
  sumstats_cache = paste0(tmpdir,'/',sha256(paste(c(p1,p2),collapse='.')),'.sumstats.rdata')
  if (args$gwas & file.exists(sumstats_cache)) load(sumstats_cache)
  if (args$gwas & ! file.exists(sumstats_cache)) {ss = sumstats(
    files = paste0(args$full, '/', c(args$p1, args$p2), '.fastGWA'), ref = args$ref,
    trait.names = trait.names, OLS = metadata$ols, N = metadata$n,
    se.logit = metadata$selogit, linprob = metadata$linprob)
    save(ss, file = sumstats_cache)
  }
  cat('\nCache files are saved at', ldsc_cache,'\n')
  if (args$gwas == 23) ss = ss %>% filter(CHR %in% c(23, 'X', '23')) else if
    (args$gwas) ss = ss %>% filter(CHR == args$gwas)
  
  #### common factor model ####
  results = paste0(args$out,'_common_mdl.txt')
  gwa = paste0(args$out,'_common.chr',args$gwas,'.fastGWA')
  if (args$common & length(trait.names_med) > 2){
    print('Compiling common factor model ...')
    if (length(p2) > 0) warning(paste0(args$p2[1],' is also included in common factor'))
    if (length(cov)> 0) warning('Covariates are also included in common factor') 
    dwls = commonfactor(ldscoutput, 'DWLS')
    out = dwls$results %>% rename(beta_unstd = 'Unstandardized_Estimate',
      se_unstd = 'Unstandardized_SE', beta = 'Standardized_Est', se = 'Standardized_SE',
      p = 'p_value') %>% mutate(q = p.adjust(p)) %>% 
      add_column(chi2 = dwls$modelfit$chisq, chi2_p = dwls$modelfit$p_chisq, 
        df = dwls$modelfit$df, CFI = dwls$modelfit$CFI, SRMR = dwls$modelfit$SRMR)
    write_tsv(out, file = results, quote = 'needed')
    
    #### common factor GWAS ####
    if (args$gwas & (!file.exists(gwa) | args$force)){
      print('Compiling common factor GWAS ...')
      if (length(args$p2) > 0) warning(paste0(args$p2[1],' is also included in common factor GWAS'))
      cgwas = commonfactorGWAS(ldscoutput, ss, cores = 16, toler = .Machine$double.xmin) %>% 
        rename(POS = 'BP', BETA = 'est', SE = 'se_c',Z = 'Z_Estimate',P = 'Pval_Estimate') %>% 
        select(-lhs, -op, -rhs, -i) %>% add_column(N = 1)
      save(cgwas, file = paste0(args$out,'_common.chr',args$gwas,'.rdata'))
      write_tsv(cgwas, gwa)
    }
  } else if (args$common & length(trait.names_med) == 2) {
    # enforce equal weight if only two factors
    print('Compiling common factor model, enforcing equal weight to two indicators ...')
    mdl = c(paste0('F1 =~ a*', trait.names_med[1],'+ a*', trait.names_med[2]),
      'F1 ~~ 1 * F1', 'a > 0.0001')
    mdl_snp = 'F1 ~ SNP'
    dwls = usermodel(ldscoutput, 'DWLS', paste0(mdl, collapse = ' \n '), 
      CFIcalc = T, imp_cov = T)
    write_model(dwls, results)
    if (args$gwas & (!file.exists(gwa) | args$force)){
      print('Compling common factor GWAS, enforcing equal weight to two indicators ...')
      mgwas = userGWAS(ldscoutput, ss, 'ML', cores = 16, smooth_check = T, sub = mdl_snp,
        model = paste0(c(mdl, mdl_snp), collapse = ' \n '), std.lv = T, toler = .Machine$double.xmin)
      write_tsv(mgwas[[1]] %>% add_column(N = 1) %>% rename(POS = 'BP', 
        BETA = 'est', Z = 'Z_Estimate',P = 'Pval_Estimate'), gwa[1])
    }
  }
  
  #### EFA and CFA ####
  results = paste0(args$out,'_cfa.txt')
  if (args$efa){
    print('Performing exploratory factor analysis ...')
    ldsccache_even = paste0(tmpdir,'/',sha256(paste(p1,collapse='.')),'.ldsc.even.rdata')
    ldsccache_odd = paste0(tmpdir,'/',sha256(paste(p1,collapse='.')),'.ldsc.odd.rdata')
    if (! file.exists(ldsccache_even)) {ldscoutput_even = ldsc(
      traits = paste0(args$input,'/',args$p1,'.sumstats'),
      sample.prev = metadata$sample_prev[1:length(p1)],
      population.prev = metadata$pop_prev[1:length(p1)],
      ld = args$ld, wld = args$ld, trait.names = p1, select = 'EVEN',
      ldsc.log = paste0(tmpdir,'/',paste(rand_bytes(4), collapse = '')))
      save(ldscoutput_even, file = ldsccache_even)} else load(ldsccache_even)
    if (! file.exists(ldsccache_odd)) {ldscoutput_odd = ldsc(
      traits = paste0(args$input,'/',args$p1,'.sumstats'),
      sample.prev = metadata$sample_prev[1:length(p1)],
      population.prev = metadata$pop_prev[1:length(p1)],
      ld = args$ld, wld = args$ld, trait.names = p1, select = 'ODD',
      ldsc.log = paste0(tmpdir,'/',paste(rand_bytes(4), collapse = '')))
      save(ldscoutput_odd, file = ldsccache_odd)} else load(ldsccache_odd)
    if (args$efa_n == -1) n_factors = paLDSC_mod(ldscoutput_even$S, ldscoutput_even$V
      ) else n_factors = args$efa_n
    s_smooth = (Matrix::nearPD(ldscoutput_even$S))$mat %>% as.matrix()
    if (n_factors == 0) n_factors = 1
    efa = factanal(covmat=s_smooth, factors = n_factors, rotation = 'promax',
      control = list(nstart = 100))
    loadings = efa$loadings
    loadings = matrix(as.numeric(loadings), attributes(loadings)$dim, 
      dimnames = attributes(loadings)$dimnames)
    colnames(loadings) = paste0('f',1:n_factors)
    loadings %>% as_tibble() %>% add_column(phenotype = p1, .before = 1) %>% 
      write_tsv(paste0(args$out, '_efa.txt'))
    
    # compile factor model formula
    print('Performing confirmatory factor analysis ...')
    mdl = paste0(p1,' ~~ var_',p1,'*',p1,' \n var_',p1,' > 0.00001')
    gwa = character(0)
    mdl_snp = character(0)
    for (i in 1:ncol(efa$loadings)) {
      eqn = character(0); used = character(0)
      for (j in 1:nrow(efa$loadings)) {
        if (abs(efa$loadings[j, i]) > args$efa_thr) eqn = c(eqn, p1[j])
      }
      if (length(eqn) == 0) eqn = p1[which.max(abs(efa$loadings[,i]))]
      eqn = paste0('F',i, ' =~ NA*', paste(eqn, collapse=' + ')); mdl = c(mdl,eqn)
      mdl_snp = c(mdl_snp, paste0('F', i,' ~ SNP'))
      gwa = c(gwa, paste0(args$out, '_efa.f',i,'.chr',args$gwas,'.fastGWA'))
    }
    dwls = usermodel(ldscoutput_odd, 'DWLS', paste0(mdl, collapse = ' \n '), 
                     CFIcalc = T, std.lv = T, imp_cov = T)
    write_model(dwls, results)
    
    #### factor GWAS ####
    if (args$gwas & (!file.exists(gwa[1]) | args$force)) {
      mgwas = userGWAS(ldscoutput, ss, 'ML', cores = 16, smooth_check = T, sub = mdl_snp,
        model = paste0(c(mdl, mdl_snp), collapse = ' \n '), std.lv = T)
      for (i in 1:length(gwa)) mgwas[[i]] %>% add_column(N = max(metadata$n)) %>% 
        rename(POS = 'BP', BETA = 'est', Z = 'Z_Estimate',P = 'Pval_Estimate') %>%
        select(-label) %>% write_tsv(gwa[i])
  }}
  
  #### exposure-outcome model ####
  results = paste0(args$out,'_causal_',length(c(med,cov)),'cov_',
    c(med,cov)[1],'.txt')
  results_med = paste0(args$out, '_',length(args$med),'med_',med[1],'.txt')
  results_adj = paste0(args$out,'_adjust_',length(cov),'cov_',cov[1],'.txt')
  results_sub = paste0(args$out, '_subtraction_mdl.txt')
  gwa = paste0(args$out,'_subtraction.chr',args$gwas,'.fastGWA')
  if (args$mdl) {
    #### direct causal model ####
    if (!file.exists(results) | args$force){
      mdl = paste0(p2[1], ' ~ NA*', paste(c(p1,cov,med), collapse = ' + '))
      dwls = usermodel(ldscoutput, model = mdl, imp_cov = T, CFIcalc = T)
      write_model(dwls, results)
    }
    
    #### adjusting for covariates ####
    if (length(cov) > 0 & (!file.exists(results_adj) | args$force)) {
      # model latent factor outcome - covars
      indep = paste0('indep_',p1)
      mdl = c(paste0('shared =~ NA*',paste(c(p1, cov), collapse = '+')), 
              paste0(indep,' =~ NA*',p1),
              'shared ~~ 1*shared',paste0(indep,'~~ 1*',indep),paste0('shared ~~ 0*',indep),
              paste0(p2, '~ NA*shared + ',indep)
              )
      zero_cov = paste0(combn(c(p1,cov),2)[1,],' ~~ 0*', combn(c(p1,cov),2)[2,])
      zero_cov = c(zero_cov, paste0(c(p1,cov),' ~~ 0*',c(p1,cov)))
      mdl = paste(c(mdl, zero_cov), collapse = '\n')
      dwls = usermodel(ldscoutput, model = mdl, imp_cov = T, CFIcalc = T)
      write_model(dwls, results_adj)
    }
    
    #### mediation model ####
    if (length(args$med) > 0 & (!file.exists(results_med) | args$force)){
      # only account for direct effect on outcome by covars, not from covars to mediators
      mdl = c(paste0(p2,'~NA*', paste(c(p1,cov,med),collapse = '+')),
        paste0(med,'~NA*',paste(p1, collapse = '+'))) %>% paste0(collapse = '\n')
      dwls = usermodel(ldscoutput, model = mdl, imp_cov = T, CFIcalc = T)
      write_model(dwls, results_med)
      if (length(cov) > 0) {
        # model with all covariates
        mdl = c(paste0(p2,'~NA*', paste(c(p1,cov,med),collapse = '+')),
                paste0(med,'~NA*',paste(c(p1,cov), collapse = '+')),
                paste0(c(p1,cov,med),' ~~ var_',c(p1,cov,med),'*',c(p1,cov,med),
                       ' \n var_',c(p1,cov,med),' > 0.00001')
        ) %>% paste0(collapse = '\n')
        dwls = usermodel(ldscoutput, model = mdl, imp_cov = T, CFIcalc = T)
        write_model(dwls, paste0(args$out,'_mediation_fullcov.txt'))
      }
    }
    
    #### subtraction model ####
    if (!file.exists(results_sub) | args$force){
      indep = paste0('indep_',p2)
      mdl = c(paste0('shared =~ NA*',paste(c(p2,cov,p1), collapse = ' + start(0.4)*')),
        paste0(indep,' =~ NA*',p2), 
        'shared ~~ 1*shared',paste0(indep,' ~~ 1*',indep),paste0('shared ~~ 0*',indep))
      zero_cov = paste0(combn(c(p1,cov,p2),2)[1,],' ~~ 0*', combn(c(p1,cov,p2),2)[2,])
      zero_cov = c(zero_cov, paste0(c(p1,cov,p2),' ~~ 0*',c(p1,cov,p2)))
      mdl = c(mdl, zero_cov) %>% paste(collapse = '\n')
      dwls = usermodel(ldscoutput, model = mdl, imp_cov = T, CFIcalc = T)
      if (is.null(dwls$results)) dwls = usermodel(ldscoutput, 'ML',
        model = mdl, imp_cov = T, CFIcalc = T)
      write_model(dwls, results_sub)
    }
    
    #### subtraction GWAS ####
    if (args$gwas & (!file.exists(gwa) | args$force)){
      indep = paste0('indep_',p2)
      mdl = c(paste0('shared =~ NA*',p2,' + start(0.2)*',
                     paste(c(p2,p1), collapse = ' + start(0.4)*')),
        paste0(indep,' =~ NA*',p2,' + start(0.2)*',p2),
        'shared ~~ 1*shared', paste0(indep,' ~~ 1*',indep),paste0('shared ~~ 0*',indep),
        'shared ~ SNP',paste0(indep,' ~ SNP'),'SNP ~~ SNP')
      zero_cov = paste0(combn(c(p1,p2),2)[1,],' ~~ 0*', combn(c(p1,p2),2)[2,])
      zero_cov = c(zero_cov, paste0(c(p1,p2),' ~~ 0*',c(p1,p2)))
      mdl = c(mdl, zero_cov) %>% paste(collapse = '\n')
      sgwas = userGWAS(ldscoutput, ss, model = mdl, cores = 16, sub = paste0(indep,' ~ SNP'))
      sgwas[[1]] %>% add_column(N = max(metadata$n)) %>% rename(POS = 'BP', BETA = 'est', 
        Z = 'Z_Estimate',P = 'Pval_Estimate') %>% select(-label) %>% write_tsv(gwa)
    }
  }
  
  #### manually input model ####
  if (!is.null(args$manual)) results = paste0(
    gsub('.mdl$','',args$manual), '.txt')
  # NB models for genome-wide model and GWAS are incompatible!
  if (!is.null(args$manual) & !args$gwas & (!file.exists(results) | args$force)) {
    mdl = read_file(args$manual); cat(mdl)
    dwls = usermodel(ldscoutput, model = mdl, imp_cov = T, CFIcalc = T)
    write_model(dwls, results)
  }
  if (!is.null(args$manual) & args$gwas) {
    mdl = read_file(args$manual); cat(mdl)
    sub = str_split(mdl,'\n')[[1]]; sub = sub[which(!is.na(str_match(sub,'SNP')))]
    sub_gwa = paste0(args$out,'_', sub %>% str_split_i('~',1) %>% gsub(' ','',.),
                     '.chr',args$gwas,'.fastGWA')
    if (!all(sapply(sub_gwa, file.exists)) | args$force){
      mgwas = userGWAS(ldscoutput, ss, model = mdl, sub = sub, cores = 16) 
      save(mgwas,file = paste0(args$out,'.usergwas.chr',args$gwas,'.rdata'))
      for (i in 1:length(sub)) mgwas[[i]] %>% add_column(N = max(metadata$n)) %>% 
        rename(POS = 'BP', BETA = 'est', Z = 'Z_Estimate',P = 'Pval_Estimate') %>%
        select(-label) %>% write_tsv(sub_gwa[i])
    }
  }
}

main(args)