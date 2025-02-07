#### Information ####
# a stand-alone script for genomic SEM analysis to analyse
# the structural contribution to the genetic correlation between functional
# graph phenotypes and autism

#### Required packages ####
library(tidyverse)
library(GenomicSEM)
library(here)

#### input files ####
func = c('../gene_corr/ldsc_sumstats/global/eff_global_0.01.sumstats.gz',
         '../gene_corr/ldsc_sumstats/global/deg_global_0.01.sumstats.gz',
         '../gene_corr/ldsc_sumstats/global/degi_global_0.01.sumstats.gz',
         '../gene_corr/ldsc_sumstats/global/degc_global_0.01.sumstats.gz',
         '../gene_corr/ldsc_sumstats/global/mpl_global_0.01.sumstats.gz',
         '../gene_corr/ldsc_sumstats/global/clu_global_0.01.sumstats.gz',
         '../gene_corr/ldsc_sumstats/global/smw_global_0.01.sumstats.gz') %>% here()
func.names = c('Global efficiency','Global degree','Ipsilateral degree',
               'Contralateral degree','Global path length','Global clustering',
               'Global small-worldness')
func.N = array(54030, dim = 7)

strc = c('../gene_corr/ldsc_sumstats/structural_for_mr/LGI_meta.sumstats.gz',
         '../gene_corr/ldsc_sumstats/structural_for_mr/Volume_meta.sumstats.gz') %>% here()
strc.names = c('Local gyrification index','Cortical volume')
strc.N = array(36663,2)

asd = here('../gene_corr/ldsc_sumstats/disorders_for_mr/asd2019.sumstats.gz')
asd.names = c('Autism')
asd.N = 46350 # NCa 18381, NCo 27969

hm3 = here('../params/ldsc_for_gsem/w_hm3.snplist')
ld = here('../toolbox/ldsc/baseline/')

