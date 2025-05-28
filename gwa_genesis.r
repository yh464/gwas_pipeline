#### Information ####
# Conducts GWAS analysis by GENESIS
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2025-05-26

#### Command line input ####
library(argparse)
parser = ArgumentParser('This script conducts GWAS by GENESIS')
parser$add_argument('pheno', help = 'Phenotype name')
parser$add_argument('-i','--in', dest = 'input', help = 'Input phenotype file, first two columns are FID and IID')
parser$add_argument('-b','--bfile', help = 'PLINK binary prefix, should have run grm_pcrelate.r')
parser$add_argument('--dcov', help = 'Discrete covariates file, first two columns are FID and IID')
parser$add_argument('--qcov', help = 'Quantitative covariates file, first two columns are FID and IID')
parser$add_argument('-o','--out', help = 'Output directory')
parser$add_argument('-f','--force', help = 'Force overwrite', default = F, action = 'store_true')
args = parser$parse_args(commandArgs(T))
args$out = normalizePath(args$out); args$bfile = normalizePath(args$bfile)
print('Input options')
print(args)

#### main execution block ####
main = function(args){
  library(tidyverse)
  library(GENESIS)
  library(SNPRelate)
  library(GWASTools)
  
  #### check progress ####
  if (! dir.exists(args$out)) dir.create(args$out)
  out_prefix = paste0(args$out, '/',args$pheno)
  if (file.exists(paste0(out_prefix, '.fastGWA')) & ! args$force) return(NULL)
  
  #### check for existence of GRM file ####
  grm = paste0(args$bfile,'.pcrelate.rdata')
  if (! file.exists(grm)){
    library(httr)
    source('https://github.com/yh464/genetics_qc/raw/refs/heads/main/grm_pcrelate.r')
    main(list(bfile = args$bfile, force = F))
  }
  load(grm)
  
  #### Process input phenotype data ####
  subj = read_table(paste0(args$bfile,'.fam'), col_names = F)
  colnames(subj) = c('FID','IID','A','B','C','D')
  subj = subj %>% select('FID','IID')
  phen = read_tsv(args$input); colnames(phen)[1:2] = c('FID','IID')
  phen = phen[,c('FID','IID',args$pheno)]
  dcov = read_tsv(args$dcov); colnames(dcov)[1:2] = c('FID','IID')
  qcov = read_tsv(args$qcov); colnames(qcov)[1:2] = c('FID','IID')
  # scanannot checks for all columns containing 'sex' to be M/F (annoying!)
  covars = c(colnames(dcov)[3:ncol(dcov)], colnames(qcov)[3:ncol(qcov)]) %>% gsub('sex','xxx',.)
  scanannot = phen %>% inner_join(subj) %>% inner_join(dcov) %>% 
    inner_join(qcov) %>% select(-FID) %>% rename(scanID = IID)
  cat('Found',nrow(scanannot),'individuals with genetic and phenotypic data')
  colnames(scanannot) = gsub('sex','xxx', colnames(scanannot))
  scanannot = ScanAnnotationDataFrame(scanannot)
  
  #### Fit null model ####
  geno = GdsGenotypeReader(paste0(args$bfile, '.gds'))
  genodata = GenotypeData(geno)
  nullmod = fitNullModel(scanannot, outcome = args$pheno, covars = covars,
                         cov.mat = pcrelategrm, family = 'gaussian')
  cat('Fit null model, time =', proc.time()[3],'\n')
  
  #### GWAS Association Test ####
  genoiter = GenotypeBlockIterator(genodata, snpBlock = 5000)
  assoc = assocTestSingle(genoiter, null.model = nullmod, 
    BPPARAM = BiocParallel::MulticoreParam(workers=16))
  h2 = varCompCI(nullmod, prop = T)
  cat('Finished association analysis, time =', proc.time()[3],'\n')
  save(nullmod, assoc, h2, file = paste0(out_prefix,'.rdata'))
  alleles = tibble(SNP = getSnpID(genodata), A1 = getAlleleA(genodata), 
                   A2 = getAlleleB(genodata))
  assoc %>% select(-MAC) %>% rename(SNP = variant.id, CHR = chr, POS = pos, 
    N = n.obs, AF1 = freq, BETA = Est, SE = Est.SE, P = Score.pval) %>% 
    filter(AF1 > 0.01 & AF1 < 0.99) %>% inner_join(alleles)
    write_tsv(paste0(out_prefix, '.fastGWA'))
  cat('Wrote GWAS output, time =', proc.time()[3],'\n')
}

main(args)