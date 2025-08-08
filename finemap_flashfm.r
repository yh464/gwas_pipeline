#### Information ####
# A wrapper for Multi-trait fine-mapping using flashfm
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2025-07-28
# Notes:  This script takes MANY summary statistics and selects SNPs in a given
#         start-stop region. Running time may be very long

#### parsing command line input ####
library(argparse)
library(here)
parser = ArgumentParser(description = 'This script runs flashfm')
# path specs
parser$add_argument('pheno', nargs = '+',
  help = 'Exposure, format <group>/<pheno>, separated by whitespace')
parser$add_argument('-i','--in', dest = 'input', help = 'input summary stats directory',
  default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/coloc/loci')
parser$add_argument('-c','--chr', help = 'chromosome', type = 'integer')
parser$add_argument('-r','--ref', 
  help = 'reference genotype file in PLINK format, auto scans directory for files split by chromosome',
  default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/bed/autosomes')
parser$add_argument('--gcov', help = 'Genetic covariance matrix')
parser$add_argument('--start', help = 'start of locus', type = 'integer')
parser$add_argument('--stop', help = 'end of locus', type = 'integer')
parser$add_argument('-o','--out', help = 'output prefix')
parser$add_argument('-f','--force', action = 'store_true', default = F, help = 'force overwrite')
args = parser$parse_args(commandArgs(TRUE))
args$pheno = sort(args$pheno); args$out = normalizePath(args$out)
if (file.exists(paste0(dirname(args$ref),'/chr',args$chr,'.bed'))) args$ref = 
  paste0(dirname(args$ref),'/chr',args$chr)
if (length(args$pheno) > 5) stop('Flashfm only supports up to 5 traits')
print('Input options')
print(args)

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
  return(list(corx = corx, snpinfo = bim))
}

#### main execution block ####
main = function(args){
  # check progress
  if (file.exists(args$out) & !args$force) return(NULL)
  
  #### read summary stats ####
  library(tidyverse)
  pheno = gsub('/','_', args$pheno)
  cache = paste0(args$input,'/',args$pheno,'_chr',args$chr,'_',args$start,'_',args$stop,'.txt') %>%
    lapply(read_tsv) %>% setNames(pheno)
  cat('Read cache files, time =', proc.time()[3])
  
  harm = cache[[1]] %>% select(SNP, A1)
  N = numeric(0)
  snps = cache[[1]]%>% select(SNP); snps = snps$SNP
  for (i in 1:length(args$pheno)){
    # intended output format: list of GWAS with columns: rsid, chromosome, position, allel1, allele2, maf, beta, se
    df = cache[[i]] %>% arrange(POS)
    N[i] = df$N %>% max() %>% as.integer()
    df = df[df$CHR == args$chr,]
    df = df[df$POS > args$start,]
    df = df[df$POS < args$stop,]
    if ('OR' %in% colnames(df)) df$BETA = log(df$OR)
    df$MAF = df$AF1; df$MAF[df$AF1 > 0.5] = 1 - df$AF1[df$AF1 > 0.5]
    df_rev = df %>% mutate(A1 = A2, BETA = -BETA)
    df = bind_rows(df, df_rev) %>% merge(harm,.) # merge with reference genome
    snps = intersect(snps, df$SNP)
    df = df %>% select(SNP, CHR, POS, A1, A2, MAF, BETA, SE) %>%
      setNames(c('rsid','chromosome','position','allele1','allele2','maf','beta','se'))
    cache[[i]] = df
  }
  for (i in 1:length(pheno)) cache[[i]] = cache[[i]] %>% filter(rsid %in% snps)
  
  #### prepare params for flashfm ####
  # reference allele freq
  ref = extract_ref(args$ref, snps)
  harm = merge(ref$snpinfo %>% select(SNP, POS, AF1), harm) %>% arrange(POS) %>% select(-POS)
  raf = harm$AF1 %>% setNames(harm$SNP)
  raf[harm$A1 != ref$snpinfo$A1] = 1 - raf[harm$A1 != ref$snpinfo$A1]
  ybar = numeric(length(args$pheno)) %>% setNames(pheno)
  N = N %>% setNames(pheno)
  if (is.null(args$gcov)) {
    covy = diag(length(pheno))
  } else covy = read.delim(args$gcov,row.names = 1); covy = covy[pheno,pheno] %>% as.matrix()
  colnames(covy) = pheno; rownames(covy) = pheno
  
  #### Run flashfm ####
  library(flashfm)
  # faulty function in flashfm package
  out = FLASHFMwithFINEMAP(cache, ref$corx, raf, ybar, N, tempdir(), 1, covy, 0.99, 1, 
    '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/finemap')
  cat('Saving temp files to', tempdir(), '\n')
  save(out, file = paste0(args$out,'.rdata'))
  
  #### parse output ####
  out.table = list()
  for (trait in pheno) {
    pp = out$mpp.pp$PPg[[trait]] %>% as.data.frame() %>% select(-null) %>% 
      setNames('pp') %>% arrange(desc(pp)) # posterior probability of credible sets, 
    trait.table = list(); total_pp = 0 # initialise for each trait
    for (i in 1:nrow(pp)) {
      tmp = tibble(pheno = trait, pp = pp$pp[i])
      snps = rownames(pp)[i] %>% str_split_1('%')
      for (snp in snps) {
        if (snp %in% names(out$snpGroups$groups.flashfm)) {
          colname = out$snpGroups$groups.flashfm[[snp]] %>% paste(collapse = '_')
        } else colname = snp
        tmp[[colname]] = 1
      }
      trait.table[[i]] = tmp
      total_pp = total_pp + pp$pp[i]
      if (total_pp > 0.95) break
    }
    out.table[[trait]] = bind_rows(trait.table)
  }
  out.table = bind_rows(out.table)
  out.table[is.na(out.table)] = 0
  write_tsv(out.table, args$out)
}

main(args)
warnings()