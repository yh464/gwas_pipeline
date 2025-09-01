#!/usr/bin/env Rscript
#### Information ####
# Creates cell-type specificity scores for each gene from scRNA-seq data using Cepo
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2025-08-28

#### Command line input ####
library(argparse)
library(here)
parser = ArgumentParser(description = 'This script runs Cepo for single-gene scoring')
parser$add_argument('-i','--in', dest = 'input', help = 'input h5ad file', required = T)
parser$add_argument('-l','--label', help = 'Columns containing cell annotations in hdf5 dataset',
  nargs = '+', default = c('supercluster_term','cluster_id','subcluster_id','roi'))
parser$add_argument('-o','--out', help = 'Output prefix, automatically appends txt suffix')
parser$add_argument('-f','--force',dest = 'force', help = 'force overwrite',
                    default = F, action = 'store_true')
args = parser$parse_args(commandArgs(TRUE))
args$input = normalizePath(args$input); args$out = normalizePath(args$out)
print(args)

main = function(args) {
  library(anndataR)
  library(Cepo)
  library(HDF5Array)
  library(rhdf5)
  # library(bit64)
  # library(Matrix)
  library(SingleCellExperiment)
  library(scater)
  library(tidyverse)
  
  # sce = read_h5ad(args$input, as = 'SingleCellExperiment', x_mapping = 'counts') %>% logNormCounts()
  h5 = H5Fopen(args$input)
  adata = read_h5ad(args$input, as = 'HDF5AnnData') # due to memory leakage, import X separately
  x = H5SparseMatrix(args$input, group = 'X') %>% as('dgCMatrix')
  sce = adata$obs %>% as.data.frame()
  rownames(x) = rownames(adata$var); colnames(x) = rownames(sce)
  gc()
  # sce = SingleCellExperiment(assays = list(logcounts = x), colData = adata$obs, rowData = adata$var)
  # sce class: row = genes, col = cells
  ref = tibble(gene = rownames(sce))
  
  scores = list()
  for (l in args$label){
    # if (! l %in% (colData(sce) %>% colnames())) {
    #   warning(paste(l, 'annotation is missing from the dataset')); next}
    scores[[l]] = (Cepo(x, sce[[l]],workers = 4))[['stats']] %>% 
      as.data.frame() %>% rename_with(~paste0('cepo.',l,'_',.x)) %>% 
      rownames_to_column('gene') %>% full_join(ref,.) %>% select(-gene)
  }
  
  scores = bind_cols(ref, scores) %>% select(where(~!all(is.na(.))))
  scores[is.na(scores)] = 0 # null = no enrichment for dropped genes
  write_tsv(scores, paste0(args$out,'.cepo.txt'))
}

main(args)