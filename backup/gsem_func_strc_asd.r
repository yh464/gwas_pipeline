#### Information ####
# a stand-alone script for genomic SEM analysis to analyse
# the structural contribution to the genetic correlation between functional
# graph phenotypes and autism

#### Required packages ####
library(tidyverse)
library(GenomicSEM)
library(here)

#### input files ####
func = c('../gene_corr/ldsc_sumstats/global/eff_global_0.01.sumstats',
         '../gene_corr/ldsc_sumstats/global/deg_global_0.01.sumstats',
         '../gene_corr/ldsc_sumstats/global/degi_global_0.01.sumstats',
         '../gene_corr/ldsc_sumstats/global/degc_global_0.01.sumstats',
         '../gene_corr/ldsc_sumstats/global/mpl_global_0.01.sumstats',
         '../gene_corr/ldsc_sumstats/global/clu_global_0.01.sumstats',
         '../gene_corr/ldsc_sumstats/global/smw_global_0.01.sumstats') %>% here()
func.raw = c('../gwa/global/eff_global_0.01.fastGWA',
             '../gwa/global/deg_global_0.01.fastGWA',
             '../gwa/global/degi_global_0.01.fastGWA',
             '../gwa/global/degc_global_0.01.fastGWA',
             '../gwa/global/mpl_global_0.01.fastGWA',
             '../gwa/global/clu_global_0.01.fastGWA',
             '../gwa/global/smw_global_0.01.fastGWA') %>% here()
func.names = c('eff_global','deg_global','degi_global',
               'degc_global','mpl_global','clu_global',
               'smw_global')
func.N = array(54030, dim = 7)
func.prev = array(NA, dim = 7)

macro = c('../gene_corr/ldsc_sumstats/structural_for_mr/LGI_meta.sumstats',
         '../gene_corr/ldsc_sumstats/structural_for_mr/Volume_meta.sumstats') %>% here()
macro.raw = c('../gwa/structural_for_mr/LGI_meta.fastGWA',
      '../gwa/structural_for_mr/Volume_meta.fastGWA') %>% here()
macro.names = c('Local gyrification index','Cortical volume')
macro.N = array(36663,2)
macro.prev = array(NA, dim = 2)

asd = here('../gene_corr/ldsc_sumstats/disorders_for_mr/asd2019.sumstats')
asd.names = c('Autism')
asd.N = 46350 # NCa 18381, NCo 27969
asd.prev = 18381/46350

hm3 = here('../params/ldsc_for_gsem/w_hm3.snplist')
ld = here('../toolbox/ldsc/baseline/')
ref = here('../params/ldsc_for_gsem/ref.1000G.txt')
out = here('../gsem/func_macro_asd/')
if (! dir.exists(out)) dir.create(out, recursive = T)

#### Common factor for functional graph phenotypes ####
func.ldsc = ldsc(traits = func, ld = ld, wld = ld, trait.names = func.names, 
                 sample.prev = func.prev, population.prev = func.prev)
func.sumstats = sumstats(files = func.raw, ref = ref, trait.names = func.names,
  se.logit = array(F, 7), OLS = array(T,7), N = func.N)
func.factor = commonfactorGWAS(covstruc = func.ldsc, SNPs = func.sumstats)

macro.ldsc = ldsc(traits = macro, ld = ld, wld = ld, trait.names = macro.names, 
                 sample.prev = macro.prev, population.prev = macro.prev)
macro.sumstats = sumstats(files = macro.raw, ref = ref, trait.names = macro.names,
  se.logit = array(F,2), OLS = array(T,2), N = func.N)
