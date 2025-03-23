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
parser$add_argument('-n','--names', help = 'Names of phenotypes, must be of same length as gwa or left blank', nargs = '*')
parser$add_argument('-s','--snp', dest = 'snp', help = 'index SNP')
parser$add_argument('-c','--chr', dest = 'chr', type = "numeric", help = 'chromosome') # these take priority over SNP
parser$add_argument('-p','--pos', dest = 'pos', type = "numeric", help = 'genomic position') # these take priority over SNP
parser$add_argument('-l', '--ld', help = 'flanking window', type = "numeric", default = 500000)
parser$add_argument('--bed',help = 'BED files to calculate LD',
#  default = '/rds/project/rds-Nl99R8pHODQ/ref/1000g')
  default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/bed')
parser$add_argument('--plink',help = 'PLINK binary',
  default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/toolbox/plink2')
parser$add_argument('-o','--out', dest = 'out', help = 'output file name')
parser$add_argument('--build', choices = c('hg19','hg38'), default = 'hg19',help = 'Human genome assembly version')
parser$add_argument('-f','--force',dest = 'force', help = 'force overwrite',
                    default = F, action = 'store_true')
args = parser$parse_args(commandArgs(TRUE))
args$out = normalizePath(args$out)
if (! endsWith(args$out,'.pdf')) args$out = paste0(args$out, '.pdf')
print('Input options')
print(args)

main = function(args) {
  if (! dir.exists(dirname(args$out))) dir.create(dirname(args$out))
  # if (file.exists(args$out) & ! args$force) return(NULL)
  
  #### extract sumstats that contain index SNP ####
  library(tidyverse)
  library(cowplot)
  
  if (! file.exists(gsub('pdf','txt',args$out)) | args$force){
    all_sumstats = tibble()
    for (f in c(1:length(args$gwa))) {
      if (is.null(args$names)){
        prefix = basename(args$gwa[f]) %>% strsplit('.', fixed = T) %>% unlist()
        prefix = prefix[-length(prefix)] %>% paste(collapse = '.') %>% 
          gsub('_0.01','',.) %>% gsub('_meta','',.)
      } else prefix = args$names[f]
      
      tmp = read.table(args$gwa[f], header=T) %>% as_tibble() %>% column_to_rownames('SNP')
      
      if (is.null(args$chr) | is.null(args$pos)) {
        chrom = tmp[args$snp,'CHR']; pos = tmp[args$snp,'POS']
        snp = args$snp
      } else {
        chrom = args$chr; pos = args$pos; 
        snp = tmp %>% filter(CHR == chrom) %>% filter(POS==pos) %>% pull(SNP)
      }
      
      tmp = tmp %>% filter(CHR == chrom)
      tmp = tmp[between(tmp$POS, pos-args$ld, pos+args$ld),]
      tmp = tmp %>% add_column(logP = -log10(tmp$P)) %>% add_column(Phenotype = prefix) 
      tmp = tmp %>% rownames_to_column('SNP')
      all_sumstats = all_sumstats %>% bind_rows(tmp)
      print(paste0('Read input file ',f))
    }
    write.table(all_sumstats, gsub('pdf','txt',args$out), sep = '\t', row.names = F)
  } else {
    all_sumstats = read.delim(gsub('pdf','txt',args$out)) %>% as_tibble()
    if (is.null(args$chr) | is.null(args$pos)) {
      prefix = all_sumstats[1,'Phenotype']
      tmp = all_sumstats %>% filter(Phenotype == prefix) %>% filter(SNP == args$snp)
      chrom = tmp[1,'CHR'] %>% as.numeric()
      pos = tmp[1,'POS'] %>% as.numeric()
      snp = args$snp
    } else {
      chrom = args$chr; pos = args$pos
      snp = all_sumstats %>% filter(CHR == chrom) %>% mutate(dist = abs(POS - args$pos)) %>%
        filter(dist == min(dist)) %>% pull(SNP)
      snp = snp[1]
    }
  }
  print('Processed input summary stats')
  
  #### Calculate R2 ####
  vcor = gsub('pdf','vcor',args$out)
  if (!file.exists(vcor)) {
    system(paste0(args$plink,' --r2-phased --bfile ', args$bed, '/chr', chrom,
                 ' --ld-window-kb ',args$ld/1000, ' --ld-snp ', snp,
                 ' --out ', gsub('.pdf','',args$out)))
    print(paste0(args$plink,' --r2-phased --bfile ', args$bed, '/chr', chrom,
                 ' --ld-window-kb ',args$ld/1000, 
                 ' --ld-window-r2 0 --ld-snp ', snp,
                 ' --out ', gsub('.pdf','',args$out)))
  }
  r2 = read.delim(vcor) %>% as_tibble() %>% select(ID_B,PHASED_R2) %>%
    bind_rows(list(ID_B = snp, PHASED_R2 = 1)%>% as.data.frame())
  colnames(r2) = c('SNP','r2')
  all_sumstats = merge(all_sumstats, r2)
  
  #### Plot locuszoom plot ####
  # change column name
  all_sumstats[['-log10(P)']] = all_sumstats$logP
  
  library(locuszoomr)
  if (args$build == 'hg19'){
    library(EnsDb.Hsapiens.v75) # hg19
    gtloc = locus(seqname = chrom, xrange = c(pos-args$ld, pos+args$ld), ens_db = 'EnsDb.Hsapiens.v75')
  } else {
    library(EnsDb.Hsapiens.v86) # hg38
    gtloc = locus(seqname = chrom, xrange = c(pos-args$ld, pos+args$ld), ens_db = 'EnsDb.Hsapiens.v86')
  }
  gt = gg_genetracks(gtloc, maxrows = 20, text_pos = 'left')
  loc = ggplot(all_sumstats, aes(POS, .data[['-log10(P)']], colour = r2)) + geom_point() +
    theme_classic() + theme(strip.background = element_blank()) + 
    scale_colour_gradient2(low = '#abdda4', mid = '#fdae61', high = '#d7191c', midpoint = 0.5) +
    xlab(paste0('Chromosome ',chrom)) +
    geom_point(aes(POS, .data[['-log10(P)']]), data = all_sumstats %>% dplyr::filter(SNP == snp),
               shape = 18, fill = 'red')+
    facet_wrap(~Phenotype, ncol = 1)
  n_pheno = length(unique(all_sumstats$Phenotype))
  fig = plot_grid(loc, gt + theme(axis.ticks.x = element_blank(), 
                                  axis.title.x = element_blank(),
                                  axis.text.x = element_blank()),
                  nrow = 2, rel_heights = c(n_pheno,1), align = 'v')
  ggsave(gsub('.pdf','.gt.pdf', args$out), plot = fig, width = 6, height = 1.75 * (1+n_pheno))
  ggsave(args$out, plot = loc, width = 6, height = 1.75 * n_pheno)
  print('Plotted locuszoom figure')
}

main(args)
warnings()