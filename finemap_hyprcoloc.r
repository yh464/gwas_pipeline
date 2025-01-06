#### Information ####
# A wrapper for Multi-trait fine-mapping using HyPrColoc
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2024-12-03
# Notes:  This script takes MANY summary statistics and selects SNPs in a given
#         start-stop region for HyPrColoc. Running time may be very long

#### reading command line input ####
library(optparse)
library(here) # for portability
optlist = list(
  # input options
  make_option(c('-i','--in'), dest = 'input', 
              help = 'input GWAS summary stats, separated by colons; requires CHR, SNP, POS, BETA/OR and SE columns'),
  make_option('--chr',dest = 'chr', type = 'integer', help = 'Chromosome'),
  make_option('--start',dest = 'start', type = 'integer', help = 'Start BP'),
  make_option('--stop',dest = 'stop', type = 'integer', help = 'Stop BP'),
  
  # output file
  make_option(c('-o','--out'),dest = 'out', help = 'output file (table)'),
  
  make_option(c('-f','--force'), dest = 'force', help = 'force overwrite', 
              default = F, action = 'store_true')
)
args = parse_args(OptionParser(option_list = optlist))
args$input = strsplit(args$input, ':')[[1]]
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
cache_betas = gsub('.txt','_betas.txt',args$out)
cache_ses = gsub('.txt','_ses.txt',args$out)

if (! file.exists(cache_betas) | ! file.exists(cache_ses) | args$force) {
  # initialise
  betas = matrix(ncol = 1, nrow = 0)
  colnames(betas) = 'SNP'
  betas = as_tibble(betas)
  ses = betas
  
  for (file in args$input){
    # intended output format: betas matrix, rows = SNP, columns = traits
    trait = paste0(basename(dirname(file)), '_',
                   basename(file) %>% gsub('.fastGWA','',.))
    print(paste0('Reading file ', file))
    df = read.table(file, header = T) %>% as_tibble()
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
    
    # clear memory
    rm(df, df1)
    
    toc = proc.time()
    print(paste0('Read ',file,', time = ',toc[3]))
  }
  write.table(betas, cache_betas, sep = '\t', quote = F)
  write.table(ses, cache_ses, sep = '\t', quote = F)
}
else {
  betas = read.table(cache_betas) %>% as_tibble()
  ses = read.table(cache_ses) %>% as_tibble()
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