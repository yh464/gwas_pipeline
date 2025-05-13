#### Information ####
# A wrapper for testing Mendelian Randomisation causality using latent
# heritable confounder (LHC-MR) and latent causal variable (LCV)
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2024-11-04
# Note:   Because the LCV package sources a broken link, the file from github
#         repo lukejoconnor/LCV is copied verbatim in here

#### The run_lcv function is copied verbatim from GitHub ####
# the 'source' command is changed
run_lcv <- function(ell,z.1,z.2,no.blocks=100,crosstrait.intercept=1,ldsc.intercept=1,
                    weights=1/pmax(1,ell),sig.threshold=.Machine$integer.max,n.1=1,n.2=1,intercept.12=0){
  mm=length(ell)
  if(length(z.1)!=mm||length(z.2)!=mm){
    stop('LD scores and summary statistics should have the same length')
  }
  source("https://github.com/lukejoconnor/LCV/raw/refs/heads/master/R/MomentFunctions.R") # original github function sourses a wrong location
  grid<- (-100:100)/100;
  # Jackknife estimates of moments
  size.blocks=floor(mm/no.blocks)
  jackknife=matrix(0,no.blocks,8)
  for(jk in 1:no.blocks){
    if(jk==1)
    {ind<-(size.blocks+1):mm}
    else if(jk==no.blocks)
    {ind <- 1:(size.blocks*(jk-1))}
    else
    {ind<-c(1:((jk-1)*size.blocks), (jk*size.blocks+1):mm)}
    jackknife[jk,] <- EstimateK4(ell[ind],z.1[ind],z.2[ind],crosstrait.intercept,ldsc.intercept,weights[ind],sig.threshold,n.1,n.2,intercept.12,8)
    if(any(is.nan(jackknife))){
      stop('NaNs produced, probably due to negative heritability estimates. Check that summary statistics and LD scores are ordered correctly.')
    }
  }
  rho.est<-mean(jackknife[,1])
  rho.err=sd(jackknife[,1])*sqrt(no.blocks+1)
  flip=sign(rho.est)
  jackknife[,2:3]<-jackknife[,2:3]-3*jackknife[,1]
  # Likelihood of each gcp value
  gcp.likelihood=grid;gcp.likelihood[]=0
  for(kk in 1:length(grid)){
    xx<-grid[kk]
    fx<-abs(jackknife[,1])^(-xx)
    numer<-jackknife[,2]/fx-fx*jackknife[,3]
    denom=pmax(1/abs(jackknife[,1]),sqrt(jackknife[,2]^2/fx^2 + jackknife[,3]^2*fx^2 ))
    pct.diff<-numer/denom
    gcp.likelihood[kk]<-dt(mean(pct.diff)/sd(pct.diff)/sqrt(no.blocks+1),no.blocks-2)
    if(xx==-1){
      pval.fullycausal.2<-pt(-flip*mean(pct.diff)/sd(pct.diff)/sqrt(no.blocks+1),no.blocks-2)
    }
    if(xx==1){
      pval.fullycausal.1<-pt(flip*mean(pct.diff)/sd(pct.diff)/sqrt(no.blocks+1),no.blocks-2)
    }
    if(xx==0){
      zscore<- flip*mean(pct.diff)/sd(pct.diff)/sqrt(no.blocks+1)
    }
  }
  pval.gcpzero.2tailed=pt(-abs(zscore),no.blocks-1)*2
  gcp.pm<-WeightedMean(grid,gcp.likelihood)
  gcp.pse<-sqrt(WeightedMean(grid^2,gcp.likelihood)-gcp.pm^2)
  h2.zscore.1<-mean(jackknife[,5])/sd(jackknife[,5])/sqrt(no.blocks+1)
  h2.zscore.2<-mean(jackknife[,6])/sd(jackknife[,6])/sqrt(no.blocks+1)
  if(h2.zscore.1<4 || h2.zscore.2<4){
    warning('Very noisy heritability estimates potentially leading to false positives')
  }
  else{
    if(h2.zscore.1<7 || h2.zscore.2<7){
      warning('Borderline noisy heritability estimates potentially leading to false positives')
    }
  }
  if(abs(rho.est/rho.err)<2){
    warning('No significantly nonzero genetic correlation, potentially leading to conservative p-values')
  }
  lcv.output<-list(zscore=zscore,pval.gcpzero.2tailed=pval.gcpzero.2tailed,gcp.pm=gcp.pm,gcp.pse=gcp.pse,rho.est=rho.est,rho.err=rho.err,
                   pval.fullycausal=c(pval.fullycausal.1,pval.fullycausal.2),h2.zscore=c(h2.zscore.1,h2.zscore.2))
  return(lcv.output)
}

#### reading command line input ####
library(optparse)
library(here) # for portability
optlist = list(
  # input options
  make_option('--g1', dest = 'gwa1', help = 'raw IDP summary stats, one file only'),
  make_option('--g2', dest = 'gwa2', help = 'raw disorder summary stats, one file only'),
  
  # LDSC scores, needs for Genomic SEM
  make_option('--ldsc', help = 'LD score file (L2), better independent from study cohorts',
              default = here('../params/ldsc_for_gsem/uk10k.l2.ldscore')), # downloaded following link from LHC website
  
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

#### infer z-score from GWAS sum stats ####
read.zscore = function(gwa){
  df = read.table(gwa, header = T)
  out = list()
  out$n = max(df$N)
  if ('OR' %in% colnames(df)) df$BETA = log(df$OR)
  if (!'Z' %in% colnames(df)) df$Z = df$BETA/df$SE # manually calculate z-score if not explicit
  out$gwa =na.omit(df[,c('SNP', 'CHR','POS','A1','A2','Z')])
  return(out)
}

#### Main LCV function ####
main = function(args) {
  library(tidyverse)
  
  #### input processing ####
  # file name operations
  prefix1 = basename(args$gwa1) %>% gsub('.?fastGWA','',.) %>% gsub('.?txt','',.)
  prefix2 = basename(args$gwa2) %>% gsub('.?fastGWA','',.) %>% gsub('.?txt','',.)
  out_prefix = paste0(args$out,'/',basename(dirname(args$gwa1)),'_',prefix1,'_',prefix2,'_mr_lcv_results.txt')
  
  # cache file specification
  tmpdir = here('../temp/mr_cache')
  if (! dir.exists(tmpdir)) dir.create(tmpdir)
  cache_file = paste0(
    tmpdir,'/',
    basename(dirname(args$gwa1)),'_',
    prefix1,'_',
    basename(dirname(args$gwa2)),'_',
    prefix2,'_lcv_cache.rdata'
  )
  
  if (! file.exists(cache_file) | args$force) {
    # read input
    gwa1 = read.zscore(args$gwa1); gwa2 = read.zscore(args$gwa2)
    n1 = gwa1$n; g1 = gwa1$gwa; n2 = gwa2$n; g2 = gwa2$gwa
    l2 = read.table(args$ldsc, sep = ',', header = T)
    if ('rs' %in% colnames(l2)) l2$SNP = l2$rs
    if ('LDSC' %in% colnames(l2)) l2$L2 = l2$LDSC
    l2 = l2[,c('SNP', 'L2')]
    
    # merge SNPs
    refsnp = l2$SNP
    refsnp = g1$SNP[g1$SNP %in% refsnp]
    refsnp = g2$SNP[g2$SNP %in% refsnp]
    refsnp = unique(refsnp)
    g1 = g1[g1$SNP %in% refsnp,]
    g1 = g1[order(g1$CHR, g1$POS),]
    g1 = g1[!duplicated(g1),]
    g2 = g2[g2$SNP %in% refsnp,]
    g2 = g2[order(g2$CHR, g2$POS),]
    g2 = g2[!duplicated(g2),]
    l2 = merge(l2, g1, by='SNP')
    l2 = l2[order(l2$CHR, l2$POS),]
    
    # harmonise alleles
    mismatch = which(g1$A1 != g2$A1, arr.ind = T)
    g2$A1[mismatch] = g1$A2[mismatch]
    g2$A2[mismatch] = g1$A1[mismatch]
    g2$Z[mismatch] = -g2$Z[mismatch]
    
    save('g1','g2','n1','n2','l2','args', file = cache_file)
  } else load(cache_file)
  
  toc = proc.time() 
  print(paste0('Finished input processing, time = ', toc[3]))
  
  if (! file.exists(out_prefix) | args$force) {
    #### run LCV ####
    res = run_lcv(ell= as.numeric(l2$L2), as.numeric(g1$Z), as.numeric(g2$Z), 
                  n.1 = as.numeric(n1), n.2 = as.numeric(n2),
                  ldsc.intercept = 0) # otherwise h2 will be underestimated and be negative
    toc = proc.time() 
    print(paste0('Finished MR-LCV, time = ', toc[3]))
    
    out = c(prefix1, prefix2, res$gcp.pm, res$gcp.pse, res$zscore, res$pval.gcpzero.2tailed,
            res$pval.fullycausal[1], res$pval.fullycausal[2], res$rho.est, res$rho.err)
    names(out) = c('exposure','outcome','gcp','se', 'z','p', 'p_fwd','p_rev', 'rho','se_rho')
    write.table(out, file = out_prefix, sep = '\t')
  }
}

main(args)
print(warnings())