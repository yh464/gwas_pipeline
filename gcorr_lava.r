#!/usr/bin/env Rscript
#### Information ####
# Flexible framework to run regional genetic correlation by LAVA
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2025-06-11

#### Command line input ####
library(argparse)
library(here)
parser = ArgumentParser(description = 'This script runs LAVA')
# path specs
parser$add_argument('-i','--in', dest = 'input', help = 'input summary stats directory',
                    default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/gwa/')
parser$add_argument('--p1', nargs = '+', required = T,
                    help = 'Exposure, format <group>/<pheno>, separated by whitespace')
parser$add_argument('--p2', nargs = '*',
                    help = 'Outcome, format <group>/<pheno>, leave blank if all pairs between p1 are needed')
parser$add_argument('--cov', nargs = '*', help = 'Covariates, format <group>/<pheno>')
parser$add_argument('--meta', nargs = '*', help = 'Metadata files')
parser$add_argument('--clump', nargs = '*', help = 'clumping outputs, in order to select for loci')
parser$add_argument('--all-loci', default = F, action = 'store_true', help = 'use all loci for regional correlation')
parser$add_argument('--overlap', help = 'sample overlap file from LDSC gcov_int')
parser$add_argument('--ref', help = 'Reference files directory', default = 
                      '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/ref/1000g_by_eth/')
parser$add_argument('--eth', help = 'ethnicity', choices = c('afr','amr','eas','eur','sas'), default = 'eur')
parser$add_argument('--all-exp', dest = 'all_exp', help = 'Multiple regression with all exposures',
                    action = 'store_true', default = F)
parser$add_argument('-o','--out', help = 'Output file name', required = T)
parser$add_argument('-f','--force', help = 'force overwrite', action = 'store_true', default = F)
args = parser$parse_args(commandArgs(TRUE))
args$p1 = sort(args$p1); args$p2 = sort(args$p2); args$cov = sort(args$cov)

#### Sanity checks to command line arguments ####
args$input = normalizePath(args$input); args$out = normalizePath(args$out)
if (length(args$meta) < length(c(args$p1,args$cov, args$p2, args$med))) warning(
  'Missing metadata, assuming traits to be continuous')
print('Input options')
print(args)

if (length(args$p2) == 0) print('Correlating all pairs of phenotypes in $p1 option') else
  if (length(args$p2) > 1) stop('Max one outcome accepted') else 
  if (!args$all_exp) print('Correlating each exposure with each outcome separately') else
  print('Correlating all exposures with each outcome in multiple regression model')
if (length(args$cov) > 0) print('All correlations will be partial, conditioned on covars')

#### Utility to parse metadata ####
parse_metadata = function(file){
  library(tidyverse)
  dat = read_tsv(file)
  colnames(dat) = colnames(dat) %>% tolower()
  if (!'group' %in% colnames(dat)) dat$group = basename(dirname(file))
  if (!'pheno' %in% colnames(dat)) dat$pheno = '*' # all phenotypes in the group
  if (!'pop_prev' %in% colnames(dat)) dat$pop_prev = NA
  if (!'nca' %in% colnames(dat)) dat$nca = NA
  if (!'nco' %in% colnames(dat)) dat$nco = NA
  dat = dat %>% select(group, pheno, pop_prev, nca, nco)
  return(dat)
}

#### Normalising LAVA outputs ####
norm_univ = function(univ) if (is.null(univ)) return(NULL) else return(univ %>% 
  rename(stat = h2.obs) %>%
  mutate(se = abs(stat/qnorm(1-p/2)),
         group1 = str_split_i(phen,':',1), pheno1 = str_split_i(phen,':',2)) %>%
  mutate(group2 = group1, pheno2 = pheno1) %>%
  select(group1, pheno1, group2, pheno2, stat, se, p) %>% add_column(type = 'h2'))
norm_biv = function(biv) if (is.null(biv)) return(NULL) else return(biv %>% 
  rename(stat = rho) %>%
  mutate(se = (rho.upper-rho.lower)/2,
         group1 = str_split_i(phen1,':',1), pheno1 = str_split_i(phen1,':',2),
         group2 = str_split_i(phen2,':',1), pheno2 = str_split_i(phen2,':',2)) %>%
  select(group1, pheno1, group2, pheno2, stat, se, p) %>% add_column(type = 'rg'))
norm_pcor = function(pcor) if (is.null(pcor)) return(NULL) else return(pcor %>%
  rename(stat = pcor) %>%
  mutate(se = (ci.upper - ci.lower)/2,
         group1 = str_split_i(phen1,':',1), pheno1 = str_split_i(phen1,':',2),
         group2 = str_split_i(phen2,':',1), pheno2 = str_split_i(phen2,':',2)) %>%
  select(group1, pheno1, group2, pheno2, stat, se, p) %>% add_column(type = 'pcor'))
norm_mreg = function(mreg) if (is.null(mreg)) return(NULL) else return(mreg %>%
  rename(state = gamma) %>%
  mutate(se = (gamma.upper-gamma.lower)/2,
         group1 = str_split_i(phen1,':',1), pheno1 = str_split_i(phen1,':',2),
         group2 = str_split_i(phen2,':',1), pheno2 = str_split_i(phen2,':',2)) %>%
  select(group1, pheno1, group2, pheno2, stat, se, p) %>% add_column(type = 'mreg'))

#### Main execution block ####
main = function(args) {
  # required packages
  library(httr)
  library(tidyverse)
  library(LAVA)
  library(openssl)
  library(progress)
  library(doParallel)
  
  tmpdir = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/temp/lava'
  if (!dir.exists(tmpdir)) dir.create(tmpdir)
  if (! dir.exists(dirname(args$out))) dir.create(dirname(args$out), recursive = T)
  if (file.exists(args$out) & !args$force) return(NULL)
  
  #### parse metadata ####
  p1 = gsub('/',':',args$p1); p2 = gsub('/',':',args$p2); cov = gsub('/',':', args$cov)
  all_metadata = lapply(args$meta, parse_metadata) %>% bind_rows()
  metadata = str_split_fixed(c(args$p1, args$cov, args$p2), '/', 2) %>% 
    as_tibble() %>% setNames(c('group','pheno')) %>% 
    add_column(pop_prev = NA, nca = NA, nco = NA)
  for (i in 1:nrow(metadata)){
    group = metadata$group[i]; pheno = metadata$pheno[i]
    match = (all_metadata$group==group) & (all_metadata$pheno %in% c(pheno,'*'))
    if (length(which(match)) > 1) stop('Conflicting metadata found!')
    if (any(match)) metadata[i,3:5] = all_metadata[match,3:5]
  }
  metadata = metadata %>% 
    mutate(filename = paste0(args$input, '/',group,'/',pheno,'.fastGWA'),
           phenotype = paste0(group,':',pheno)) %>%
    rename(cases = nca, controls = nco, prevalence = pop_prev) %>%
    select(phenotype, cases, controls, prevalence, filename)
  infofile = paste0(rand_bytes(8),collapse = '') %>% paste0(tempdir(), .,'.info.txt')
  temp_files = infofile
  write_tsv(metadata, file = infofile)
  
  #### read loci ####
  blocks = list()
  source('https://github.com/cadeleeuw/lava-partitioning/raw/refs/heads/main/ldblock.r')
  for (chrom in 1:22){
    blocks[[chrom]] = load.breaks(paste0(args$ref,'/',args$eth,'/chr',chrom)) %>% 
      make.blocks() %>% as_tibble() %>% add_column(CHR = chrom, .before = 1)
  }
  blocks = bind_rows(blocks) %>% rename(START = start, STOP = stop) %>% add_column(sig = F)
  
  # filter loci based on significance
  if (length(args$clump) == 0) {
    args$all_loci = T
    warning('Missing clump file, correlating all loci')
  }
  if (! args$all_loci){
    clumps = lapply(args$clump, read_tsv) %>% bind_rows() %>% select(CHR, BP) %>% 
      arrange(CHR, BP)
    nclump = nrow(clumps); clump_id = 1; clump = clumps[1,]
    for (i in 1:nrow(blocks)){
      block = blocks[i,]
      while ((block$CHR == clump$CHR) & (block$START <= clump$BP) & (block$STOP >= clump$BP)){
        blocks$sig[i] = T
        clump_id = clump_id + 1
        if (clump_id > nclump) break
        clump = clumps[clump_id,]
      }
      if (clump_id > nclump) break
    }
    blocks = blocks %>% filter(sig)
  }
  loci = blocks %>% select(CHR,START,STOP) %>% add_column(LOC=1:nrow(blocks), .before = 1)
  cat('\n\nFound',nrow(loci),'loci for analysis\n\n')
  
  #### read input sumstats ####
  cache = paste0(tmpdir,'/',sha256(paste(p1,p2,cov, collapse = '_')),'.lava.rdata')
  if (!file.exists(cache) | args$force){
    input = process.input(
      input.info.file = infofile, sample.overlap.file = args$overlap,
      ref.prefix = paste0(args$ref,'/',args$eth, '/autosomes'),
    )
    save(input, file = cache)
  } else load(cache)
  
  #### run LAVA ####
  # registerDoParallel(makeCluster(4))
  pb = progress_bar$new(total = nrow(loci), format = '[:bar] :percent :elapsedfull / eta :eta')
  pb$tick(0)
  stats_files = character(0)
  error_loci = list()
  # stats = foreach(i=1:nrow(loci)) %dopar% {
  for (i in 1:nrow(loci)) {
    tmpfile = paste0(tmpdir,'/',sha256(paste(p1,p2,cov, collapse = '_')),'_locus_',i,'.txt')
    temp_files[i+1]=tmpfile
    if (file.exists(tmpfile) & !args$force) {stats_files[i] = tmpfile; pb$tick(); next}
    locus = tryCatch(process.locus(loci[i,], input), error = function(e) NULL)
    if (is.null(locus)) {error_loci[[i]] = loci[i,]; next}
    cat('Processed locus, Time =',proc.time()[3], '\n')
    locus_stats = list()
    
    if (length(p2) == 0) {
      # all pairs of bivariate correlations
      locus_stats$univ = run.univ(locus, phenos = p1) %>% norm_univ()
      locus_stats$biv = run.bivar(locus, phenos = p1) %>% norm_biv()
      if (length(cov) > 0) {
        # partial correlations
        pcor = list()
        pairwise = expand_grid(phen1 = p1, phen2 = p1) %>% as.matrix()
        for (j in 1:nrow(pairwise)){
          pcor[[j]] = run.multireg(locus, target = pairwise[j,], phenos = cov)}
        locus_stats$pcor = bind_rows(pcor) %>% norm_pcor()
      }
    } else if (! args$all_exp) {
      # pairwise correlation with all outcome phenotypes
      locus_stats$univ = run.univ(locus, phenos = c(p1, p2)) %>% norm_univ()
      locus_stats$biv = run.bivar(locus, phenos = c(p1, p2), target = p2) %>% norm_biv()
      if (length(cov) > 0){
        pcor = list()
        for (j in 1:length(p1)){
          pcor[[j]] = run.multireg(locus, target = c(p1[j],p2), phenos = cov)}
        locus_stats$pcor = bind_rows(pcor) %>% norm_pcor()
      }
    } else {
      locus_stats$univ = run.univ(locus, phenos = c(p1,p2)) %>% norm_univ()
      locus_stats$mreg = run.multireg(locus, target = p2, phenos = c(p1, cov, p2)) %>% norm_mreg()
    }
    locus_stats = bind_rows(locus_stats) %>% add_column(
      CHR = locus$chr, START = locus$start, STOP = locus$stop, nsnp = locus$n.snps
    )
    if (! is.null(locus_stats)){
      stats_files[[i]] = tmpfile
      write_tsv(locus_stats, tmpfile)
    } else error_loci[[i]]= loci[i,]
    cat('\n',i,'/',nrow(loci), 'time =',proc.time()[3],'\n')
    pb$tick()
  }
  stats = lapply(stats_files %>% na.omit(), read_tsv, col_types = 
    cols(CHR = 'i', START='d', STOP = 'd', nsnp = 'i')) %>% bind_rows() %>% 
    mutate(p_bonferroni = p * nrow(loci)) %>% 
    group_by(type) %>% mutate(q = p.adjust(p, 'BH'))
  write_tsv(stats, args$out)
  
  if (length(error_loci) > 0) write_tsv(bind_rows(error_loci), 
    gsub('.txt','',args$out) %>% paste0('.err'))
  # clean up
  # for (file in temp_files) file.remove(file)
}

main(args)
