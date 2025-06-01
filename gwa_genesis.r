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
parser$add_argument('--n_threads', help = 'Number of threads', type = 'integer', default = 8)
# because GENESIS is used for small populations, use UKB as a reference for MAF filtering
parser$add_argument('--snp', help = 'SNP information table of reference population, used for MAF filter',
  default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/params/ukb_snp_info.txt')
parser$add_argument('--extract', help = 'List of SNPs to extract')
parser$add_argument('--maf', help = 'MAF filter', type = 'numeric', default = 0.01)
parser$add_argument('-o','--out', help = 'Output directory')
parser$add_argument('-f','--force', help = 'Force overwrite', default = F, action = 'store_true')
args = parser$parse_args(commandArgs(T))
args$out = normalizePath(args$out); args$bfile = normalizePath(args$bfile)
print('Input options')
print(args)

#### main execution block ####
main = function(args){
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
  library(tidyverse)
  library(GENESIS)
  library(SNPRelate)
  library(GWASTools)
  subj = read_table(paste0(args$bfile,'.fam'), col_names = F)
  colnames(subj) = c('FID','IID','A','B','C','D')
  subj = subj %>% select('FID','IID')
  phen = read_tsv(args$input); colnames(phen)[1:2] = c('FID','IID')
  phen = phen[,c('FID','IID',args$pheno)]
  colnames(phen) = c('FID','IID','phen')
  dcov = read_tsv(args$dcov); colnames(dcov)[1:2] = c('FID','IID')
  qcov = read_tsv(args$qcov); colnames(qcov)[1:2] = c('FID','IID')
  # scanannot checks for all columns containing 'sex' to be M/F (annoying!)
  covars = c(colnames(dcov)[3:ncol(dcov)], colnames(qcov)[3:ncol(qcov)]) %>% gsub('sex','xxx',.)
  scanannot = phen %>% inner_join(subj) %>% inner_join(dcov) %>% 
    inner_join(qcov) %>% select(-FID) %>% rename(scanID = IID)
  cat('Found',nrow(scanannot),'individuals with genetic and phenotypic data\n')
  colnames(scanannot) = gsub('sex','xxx', colnames(scanannot))
  scanannot = ScanAnnotationDataFrame(scanannot)
  
  #### Fit null model ####
  geno = GdsGenotypeReader(paste0(args$bfile, '.gds'))
  genodata = GenotypeData(geno)
  nullmod = fitNullModel(scanannot, outcome = 'phen', covars = covars,
                         cov.mat = pcrelategrm, family = 'gaussian')
  cat('Fit null model, time =', proc.time()[3],'\n')
  
  #### Select SNPs to include ####
  if (! is.null(args$extract)) {snps = tibble(SNP = read_lines(args$extract, 
    skip_empty_rows = T)) %>% dplyr::intersect(tibble(SNP = getSnpID(genodata)))
    snps = snps$SNP
  } else if (file.exists(args$snp) & args$maf > 0){
    snps = read_tsv(args$snp) %>% filter((AF1 > args$maf) & (AF1 < 1-args$maf)) %>%
      select(SNP) %>% dplyr::intersect(tibble(SNP = getSnpID(genodata)))
    snps = snps$SNP
    cat('Including',length(snps),'SNPs for analysis\n')
    cat('Filtered MAF, time =', proc.time()[3],'\n')
  } else snps = NULL
  
  #### GWAS Association Test ####
  genoiter = GenotypeBlockIterator(genodata, snpBlock = 5000, snpInclude = snps)
  assoc = assocTestSingle(genoiter, null.model = nullmod, 
    BPPARAM = BiocParallel::MulticoreParam(workers=args$n_threads))
  h2 = varCompCI(nullmod, prop = T)
  cat('Finished association analysis, time =', proc.time()[3],'\n')
  
  alleles = tibble(SNP = getSnpID(genodata), A1 = getAlleleA(genodata), 
                   A2 = getAlleleB(genodata))
  assoc %>% as_tibble() %>% select(-MAC) %>% rename(SNP = variant.id, CHR = chr, POS = pos, 
    N = n.obs, AF1 = freq, BETA = Est, SE = Est.SE, P = Score.pval) %>% inner_join(alleles) %>%
    mutate(N = N/2) %>% write_tsv(paste0(out_prefix, '.fastGWA'))
  cat('Wrote GWAS output, time =', proc.time()[3],'\n')
}

main(args)