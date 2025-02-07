require(data.table)
require(dplyr)
require(tidyverse)
require(argparse)
require(GenomicSEM)
require(Matrix)
require(stats)

parser <- ArgumentParser(description = 'Prepares data for Genomic SEM pipeline')
parser$add_argument('--ldsc', dest = 'ldsc', help = 'dir of LDSC sumstats',
                    default = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/gene_corr/ldsc_sumstats/')
parser$add_argument('--list', dest = 'list', help =
                      'List of fastGWA files to parse')
parser$add_argument('--file', dest = 'file', nargs = '*', help =
                      'Input fastGWA file')
parser$add_argument('--ref', dest = 'ref', help = 'reference ldsc files',
                    default = '/rds/user/yh464/hpc-work/ldsc/baseline/')
parser$add_argument('--fref', dest = 'fref', help = 'reference MAF file',
                    default = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/params/g1000_eur_0.01.txt')
parser$add_argument('-o','--out', dest = 'out', help = 'output directory',
                    default = '/rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/multivar-gwa/')
parser$add_argument('-p','--prefix', dest = 'prefix', help = 'name of the output file',
                    required = T)
parser$add_argument('-f', '--force', dest = 'force', action = 'store_true',
                    default = F, help = 'force overwrite')
args <- parser$parse_args()

# directory operations
if (!dir.exists(args$out)) {dir.create(args$out)}
setwd(args$out)

# parse input list
flist <- readLines(args$list)
flist <- c(args$file, flist)
prefix <- NULL
munged <- NULL

for (x in flist) {
  if (!grep('.fastGWA', x, ignore.case = T)) {
    flist <- flist[flist!=x]
    next}
  tmp <- strsplit(x, '/')
  tmp <- tmp[[1]]
  tmp <- tmp[length(tmp)]                                                       # remove the directories
  tmp <- gsub('.fastGWA', '', tmp)
  prefix <- c(prefix, tmp)
  tmp2 <- sprintf('%s%s.sumstats',args$ldsc,tmp)
  munged <- c(munged, tmp2)
}

cat(' ')
print(prefix)
cat(' ')

# NO NEED TO MUNGE USING THE SAME LDSC FUNCTION!
# munge(files = flist, hm3 = args$ref, prefix = prefix,
#       info.filter = 0.9, maf.filter = 0.1)

# input parameters
sprev <- rep(NA, times = length(prefix))
pprev <- sprev
se.logit <- rep (F, times = length(prefix))
ols <- rep(T, times = length(prefix))
linprob <- se.logit

# file names
ldsc_fname <- sprintf('%s%s_LDSCoutput.RData', args$out, args$prefix)
mvss_fname <- sprintf('%s%s_mvsumstats.RData', args$out, args$prefix)

if (!file.exists(ldsc_fname) | args$force) {
# multivariate LDSC
LDSCoutput <- ldsc(
  traits = munged,
  trait.names = prefix,
  sample.prev = sprev, population.prev = pprev,
  ld = args$ref, wld = args$ref,
)
save(LDSCoutput, file = ldsc_fname)
}

if (!file.exists(mvss_fname) | args$force) {
# prepare sumstats
mvsumstats <- sumstats(
  files = flist,
  ref = args$fref,
  trait.names = prefix,
  se.logit = se.logit, OLS = ols, linprob = linprob,
  betas = ols, N = sprev,
  maf.filter = .01,
  keep.indel = T,
  parallel = T,
  cores = 4
)
save(mvsumstats, file = mvss_fname)
}
