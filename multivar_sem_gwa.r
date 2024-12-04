require(data.table)
require(dplyr)
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

# file names
ldsc_fname <- sprintf('%s%s_LDSCoutput.RData', args$out, args$prefix)
mvss_fname <- sprintf('%s%s_mvsumstats.RData', args$out, args$prefix)
out_fname <- sprintf('%s%s_mvgwa.txt', args$out, args$prefix)
out_rdata <- sprintf('%s%s_mvgwa.RData', args$out, args$prefix)

# step 4 - combine sumstats and LDSC output
load(ldsc_fname)
load(mvss_fname)
mvgwa = commonfactorGWAS(
  covstruc = LDSCoutput,
  SNPs = mvsumstats,
  estimation = 'DWLS',
  cores = 32,
  toler = F,
  parallel = T,
)
write.csv(mvgwa, file = out_fname, sep = '\t', row.names = F)
save(mvgwa, file = out_rdata)
