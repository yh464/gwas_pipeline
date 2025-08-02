#### Information ####
# A wrapper for Multi-trait fine-mapping using flashfm
# Author: Yuankai He (yh464@cam.ac.uk)
# Date:   2025-07-28
# Notes:  This script takes MANY summary statistics and selects SNPs in a given
#         start-stop region. Running time may be very long

#### parsing command line input ####
library(argparse)
library(here)
parser = ArgumentParser(description = 'This script runs LAVA')
# path specs
parser$add_argument('pheno', nargs = '+', required = T,
  help = 'Exposure, format <group>/<pheno>, separated by whitespace')
parser$add_argument('-i','--in', dest = 'input', help = 'input summary stats directory',
  default = '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/gwa/')
parser$add_argument('-c','--chr', help = 'chromosome', type = numeric)
parser$add_argument('--start', help = 'start of locus', type = numeric)
parser$add_argument('--stop', help = 'end of locus', type = numeric)
parser$add_argument('-o','--out', help = 'output prefix')
parser$add_argument('-f','--force', action = 'store_true', default = F, help = 'force overwrite')
args = parser$parse_args(commandArgs(TRUE))
args$pheno = sort(args$pheno); args$out = normalizePath(args$out)

print('Input options')
print(args)

#### extract summary stats ####
