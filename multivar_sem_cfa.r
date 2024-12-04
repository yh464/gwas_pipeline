require(data.table)
require(dplyr)
require(argparse)
require(GenomicSEM)
require(Matrix)
require(stats)

parser <- ArgumentParser(description = 'Runs genomic EFA and CFA')
parser$add_argument('--list', dest = 'list', help =
                      'List of fastGWA files to parse')
parser$add_argument('--file', dest = 'file', nargs = '*', help =
                      'Input fastGWA file')
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

for (x in flist) {
  if (!grep('.fastGWA', x, ignore.case = T)) {
    flist <- flist[flist!=x]
    next}
  tmp <- strsplit(x, '/')
  tmp <- tmp[[1]]
  tmp <- tmp[length(tmp)]                                                       # remove the directories
  tmp <- gsub('.fastGWA', '', tmp)
  prefix <- c(prefix, tmp)
}

cat(' ')
print(prefix)
cat(' ')

# file names
ldsc_fname <- sprintf('%s%s_LDSCoutput.RData', args$out, args$prefix)
mvss_fname <- sprintf('%s%s_mvsumstats.RData', args$out, args$prefix)
log_fname <- sprintf('%s%s_CFA.log', args$out, args$prefix)
sink(log_fname)

load(ldsc_fname)
ssmooth  <- as.matrix((nearPD(LDSCoutput$S, corr = F))$mat)
print(ssmooth)

fail_flag <- F

for (x in 1:nrow(ssmooth)) {
  if (x >= nrow(ssmooth) * 3/4) break
  test <- try(
    factanal(covmat = ssmooth, factors = x),
    silent = T
  )
  print(test)
  if (class(test) == 'try-error') {
    fail_flag <- T
    break
  }
  efa <- factanal(covmat = ssmooth, factors = x)
  # specify efa model
  mdl <- ''
  for (i in 1:x) {
    mdl <- paste(mdl, 'F', i, ' =~ NA*', sep = '')
    tmp <- matrix(efa$loadings)
    print(tmp)
    print(i)
    for (j in 1:length(prefix)) {                                               # extract sig. factors
      if (abs(tmp[j,i]) > 0.1) {                                       # sig. loading
        if (substr(mdl, length(mdl), length(mdl)) != '*') {
          mdl <- paste(mdl, '+', prefix[j])
        }
        else {
          mdl <- paste(mdl, prefix[j], sep = '')
        }
      }
    }
    mdl <- paste(mdl, '\n')                                                     # new line after each eqn
  }
  if (x > 1) {
    for (i in 1:x) {
      for (j in 1:i-1) {
        mdl <- paste(mdl, 'F',i, '~~F', j, '\n', sep = '')
      }
    }
  }
  print(mdl)
  cfa <- usermodel(LDSCoutput, estimation = 'DWLS', model = mdl,
                   CFIcalc = T, std.lv = T, imp_cov = F)

  print(cfa)
}

# if the matrix is singular, use a one-factor model with all traits involved
if (fail_flag) {
  mdl <- paste(prefix,collapse = ' + ')
  mdl <- paste('F1 =~ NA*', mdl, sep = '')
  print(mdl)
  cfa <- usermodel(LDSCoutput, estimation = 'DWLS', model = mdl,
                   CFIcalc = T, std.lv = T, imp_cov = F)

  print(cfa)
}

sink() # redirects R out put to stdout
