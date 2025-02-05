#!/usr/bin/env Rscript
#### Information ####
# Batch plots locuszoom plots for specified input GWAS summary stats and one single SNP
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2025-01-31

#### Parsing cmdline input ####
library(argparse)
library(here)
parser = ArgumentParser(description = 'This programme compiles LocusZoom plots for any number of input GWAS summary stats')
parser$add_argument('gwa', help = 'GWAS summary stats', nargs = '*')
parser$add_argument('-s','--snp', dest = 'snp', help = 'index SNP')
parser$add_argument('-c','--chr', dest = 'chr', type = "numeric", help = 'chromosome') # these take priority over SNP
parser$add_argument('-p','--pos', dest = 'pos', type = "numeric", help = 'genomic position') # these take priority over SNP
parser$add_argument('-l', '--ld', help = 'flanking window', type = "numeric", default = 1000000)
parser$add_argument('-o','--out', dest = 'out', help = 'output file name')
parser$add_argument('--build', choices = c('hg19','hg38'), default = 'hg19',help = 'Human genome assembly version')
parser$add_argument('-f','--force',dest = 'force', help = 'force overwrite',
                    default = F, action = 'store_true')
args = parser$parse_args(commandArgs(TRUE))
print('Input options')
print(args)

main = function(args) {
  if (! dir.exists(dirname(args$out))) dir.create(dirname(args$out))
  # if (file.exists(args$out) & ! args$force) return(NULL)
  
  # extract sumstats that contain index SNP
  library(tidyverse)
  library(cowplot)
  
  if (! file.exists(paste0(args$out, '.txt')) | args$force){
    all_sumstats = tibble()
    for (f in args$gwa) {
      prefix = basename(f) %>% strsplit('.', fixed = T) %>% unlist()
      prefix = prefix[-length(prefix)] %>% paste(collapse = '.') %>% 
        gsub('_0.01','',.) %>% gsub('_meta','',.)
      tmp = read.table(f, header=T) %>% as_tibble() %>% column_to_rownames('SNP')
      
      if (is.null(args$chr) | is.null(args$pos)) {chrom = tmp[args$snp,'CHR']; pos = tmp[args$snp,'POS']
      } else {chrom = args$chr; pos = args$pos}
      
      tmp = tmp %>% filter(CHR == chrom) %>% filter(between(POS, pos-args$ld, pos+args$ld)) 
      tmp = tmp %>% add_column(logP = -log10(tmp$P)) %>% add_column(Phenotype = prefix) 
      tmp = tmp %>% rownames_to_column('SNP')
      all_sumstats = all_sumstats %>% bind_rows(tmp)
      print(paste0('Read input file ',f))
    }
    write.table(all_sumstats, paste0(args$out, '.txt'), sep = '\t', row.names = F)
  } else {
    all_sumstats = read.table(paste0(args$out, '.txt'), header=T) %>% as_tibble()
    if (is.null(args$chr) | is.null(args$pos)) {
      prefix = all_sumstats[1,'Phenotype']
      tmp = all_sumstats %>% filter(Phenotype == prefix) %>% filter(SNP == args$snp)
      chrom = tmp[1,'CHR'] %>% as.numeric()
      pos = tmp[1,'POS'] %>% as.numeric()
    } else {chrom = args$chr; pos = args$pos}
  }
  print('Processed input summary stats')
  
  library(locuszoomr)
  if (args$build == 'hg19'){
    library(EnsDb.Hsapiens.v75) # hg19
    gtloc = locus(seqname = chrom, xrange = c(pos-args$ld, pos+args$ld), ens_db = 'EnsDb.Hsapiens.v75')
  } else {
    library(EnsDb.Hsapiens.v86) # hg38
    gtloc = locus(seqname = chrom, xrange = c(pos-args$ld, pos+args$ld), ens_db = 'EnsDb.Hsapiens.v86')
  }
  gt = gg_genetracks(gtloc, maxrows = 20, text_pos = 'left')
  loc = ggplot(all_sumstats, aes(POS, logP, colour = Phenotype)) + geom_point() +
    theme_light() + 
    theme(axis.ticks.x = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  fig = plot_grid(loc, gt, nrow = 2, rel_heights = c(1,1), align = 'v')
  ggsave(args$out, plot = fig, width = 12, height = 7)
  print('Plotted locuszoom figure')
}

main(args)
warnings()