#### Information ####
# A wrapper for Multi-trait fine-mapping using mvSuSiE
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2025-12-14
# Notes:  This script takes MANY summary statistics and selects SNPs in a given
#         start-stop region. Running time may be very long

#### parsing command line input ####
library(argparse)
library(here)
parser = ArgumentParser(description = 'This script runs mvSuSiE')
# path specs
parser$add_argument('pheno', nargs = '+',
                    help = 'Exposure, format <group>/<pheno>, separated by whitespace')
parser$add_argument('-i','--in', dest = 'input', help = 'input summary stats directory',
                    default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/coloc/loci')
parser$add_argument('-r','--ref', 
                    help = 'reference genotype file in PLINK format, auto scans directory for files split by chromosome',
                    default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/bed/chr%chr%')
parser$add_argument('--gcov', help = 'Genetic covariance matrix')
parser$add_argument('-c','--chr', help = 'chromosome', type = 'integer')
parser$add_argument('--start', help = 'start of locus', type = 'integer')
parser$add_argument('--stop', help = 'end of locus', type = 'integer')
parser$add_argument('-o','--out', help = 'output prefix')
parser$add_argument('-f','--force', action = 'store_true', default = F, help = 'force overwrite')
args = parser$parse_args(commandArgs(TRUE))
args$pheno = sort(args$pheno); args$out = normalizePath(args$out)
if (file.exists(paste0(dirname(args$ref),'/chr',args$chr,'.bed'))) args$ref = 
  paste0(dirname(args$ref),'/chr',args$chr)
print('Input options')
print(args)

#### read all sumstats for a specific locus ####
read_sumstats = function(args){
  library(tidyverse)
  pheno = gsub('/','_', args$pheno)
  cache = paste0(args$input,'/',args$pheno,'_chr',args$chr,'_',args$start,'_',args$stop,'.txt') %>%
    lapply(read_tsv) %>% setNames(pheno)
  snps = cache %>% lapply(select, SNP) %>% Reduce(intersect, .)
  snps = snps$SNP
  harm = cache[[1]] %>% select(SNP, A1)
  N = numeric(0)
  
  for (i in 1:length(args$pheno)){
    # intended output format: a harmonised N_snp * N_pheno matrix of beta and SE
    df = cache[[i]] %>% arrange(POS)
    N[i] = df$N %>% max() %>% as.integer()
    df = df[df$CHR == args$chr,]
    df = df[df$POS > args$start,]
    df = df[df$POS < args$stop,]
    if ('OR' %in% colnames(df)) df$BETA = log(df$OR)
    df_rev = df %>% mutate(A1 = A2, BETA = -BETA)
    df = bind_rows(df, df_rev) %>% merge(harm,.) # merge with reference genome
    df = df %>% filter(SNP %in% snps) %>% select(SNP, A1, A2, BETA, SE)
    cache[[i]] = df
  }
  snps = cache[[1]]$SNP # force SNPs to be in the correct order
  
  out_beta = cache %>% lapply(select, BETA) %>% bind_cols() %>% as.matrix()
  colnames(out_beta) = pheno; rownames(out_beta) = snps
  out_se = cache %>% lapply(select, SE) %>% bind_cols() %>% as.matrix()
  colnames(out_se) = pheno; rownames(out_se) = snps
  
  cat('Read cache files, time =', proc.time()[3], '\n')
  return(list(bhat = out_beta, shat = out_se, N = as.integer(max(N)), snps = snps))
}

#### extract reference genotypes ####
extract_ref = function(prefix, snp.list){
  library(tidyverse)
  library(snpStats)
  # first extract the required SNPs
  bim = read_delim(paste0(prefix, '.bim'), col_names = c('CHR','SNP','cM','POS','A1','A2')) %>% 
    filter(SNP %in% snp.list) %>% select(-CHR, -cM)
  bed = read.plink(prefix, select.snps = bim$SNP)
  corx = ld(bed$genotype, depth = nrow(bim), stats = 'R', symmetric = T) %>% as.matrix()
  # snp = rownames(corx) %>% gsub(':','_',.)
  # snp[any(startsWith(snp, as.character(1:9)))] = paste0('chr', snp[any(startsWith(snp, as.character(1:9)))])
  diag(corx) = 1
  snp.stats = col.summary(bed$genotype)
  raf = 1 - snp.stats$RAF
  bim$AF1 = raf
  
  cat('Estimated LD matrix, time =', proc.time()[3], '\n')
  return(list(ld = corx, snpinfo = bim))
}

#### main execution block ####
main = function(args){
  library(tidyverse)
  library(mvsusieR)
  
  sumstats = read_sumstats(args)
  ref = extract_ref(args$ref %>% gsub('%chr%', args$chr, .), sumstats$snps)
  # merge SNPs from PLINK ref and sumstats
  sumstats$snps = sumstats$snps[sumstats$snps %in% ref$snpinfo$SNP]
  sumstats$bhat = sumstats$bhat[sumstats$snps,]
  sumstats$shat = sumstats$shat[sumstats$snps,]
  
  res = mvsusie_rss(R = ref$ld, N = sumstats$N, Bhat = sumstats$bhat, Shat = sumstats$shat)
  print(res)
}

main(args)
warnings()