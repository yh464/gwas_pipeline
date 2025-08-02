#### Information ####
# A wrapper for Multi-trait fine-mapping using HyPrColoc
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2024-12-03
# Notes:  This script takes MANY summary statistics and selects SNPs in a given
#         start-stop region for HyPrColoc. Running time may be very long

#### reading command line input ####
library(argparse)
library(here) # for portability
parser = ArgumentParser(description = 'This script runs HyPrColoc based on extracted loci')
parser$add_argument('pheno', nargs = '*', help = 'Phenotypes in format <group>/<pheno> separated by whitespace')
parser$add_argument('-i','--in', dest = 'input', help = 'Input directory containing extracted loci',
  default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/coloc/loci')
parser$add_argument('--chr', type = 'integer', help = 'Chromosome, 23 for X, 24 for Y, 25 for XY, 26 for MT')
parser$add_argument('--start', type = 'integer', help = 'Start BP')
parser$add_argument('--stop', type = 'integer', help = 'Stop BP')
parser$add_argument('-o','--out', help = 'Output file .txt')
parser$add_argument('-f','--force', help = 'Force overwrite', default = F, action = 'store_true')
args = parser$parse_args(commandArgs(TRUE))
args$input = normalizePath(args$input); args$out = normalizePath(args$out)
print('Input options')
print(args)


main = function(args){

#### Required packages ####
if (! require(tidyverse)){
  install.packages('tidyverse', repos = "https://cloud.r-project.org")
  library(tidyverse)
}
library(hyprcoloc)

#### Read all input files and extract relevant SNPs ####
cache = paste0(args$input,'/',args$pheno,'_chr',args$chr,'_',args$start,'_',args$stop,'.txt') %>%
  lapply(read_tsv)
cat('Read cache files, time =', proc.time()[3])

# initialise
betas = matrix(ncol = 1, nrow = 0)
colnames(betas) = 'SNP'
betas = as_tibble(betas)
ses = betas

for (i in 1:length(args$pheno)){
  # intended output format: betas matrix, rows = SNP, columns = traits
  df = cache[[i]]; trait = args$pheno[i]
  df = df[df$CHR == args$chr,]
  df = df[df$POS > args$start,]
  df = df[df$POS < args$stop,]
  if ('OR' %in% colnames(df)) df$BETA = log(df$OR)
  df = df[,c('CHR','SNP','POS','BETA','SE')]
  
  # outer merge to identify problematic files
  df1 = df[,c('SNP','BETA')]
  colnames(df1) = c('SNP',trait)
  betas = merge(betas, df1, all = T)
  
  df1 = df[,c('SNP','SE')]
  colnames(df1) = c('SNP',trait)
  ses = merge(ses, df1, all = T)
}


#### Run Hyprcoloc ####
ses[ses==0] = 1e-8
betas = betas %>% drop_na() %>% column_to_rownames('SNP')
ses = ses %>% drop_na() %>% column_to_rownames('SNP')

toc = proc.time()
print(paste0('Finished input processing, time = ',toc[3]))
print(nrow(ses))
res = hyprcoloc(as.matrix(betas), as.matrix(ses), 
  trait.names = colnames(betas), snp.id = rownames(betas))
write.table(res$results, args$out, sep = '\t')
}

main(args)
print(warnings())