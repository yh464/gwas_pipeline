#!/usr/bin/env Rscript
#### Information ####
# Normalises text in R data frames
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2025-03-04

library(tidyverse)
dict = normalizePath('../path/dict.txt')
dict = read.delim(dict)

normalise_df = function(df) {
  tmp = df %>% select(!SNP)
  for (i in c(1:nrow(dict))) {
    before = dict[i,'before']; after = dict[i,'after']
    
    if (startsWith(before, '_')) {
      tmp = gsub(before, after, tmp)
      next
    } else {
      tmp = gsub(paste0('_', before, '_'), paste0('_',after,'_'), tmp)
      tmp = gsub(paste0('^', before, '_'), paste0(after, '_'), tmp)
      tmp = gsub(paste0('_', before, '$'), paste0('_',after), tmp)
      tmp = gsub(paste0('^', before, '$'), after, tmp)
    }
  }
  colnames(tmp) = colnames(df)
  rownames(tmp) = rownames(df)
  tmp = as_tibble(tmp)
  if ('SNP' %in% colnames(df)) tmp = bind_cols(tmp, df$SNP)  
  return(tmp)
}