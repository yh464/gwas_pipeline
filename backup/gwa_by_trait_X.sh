#!/bin/bash
#!/bin/bash
$1 --fastGWA-mlm --mbfile $2 --grm-sparse $3 --pheno $4 --mpheno $5 --qcovar $6 --covar $7 --out $8 ${@:10}
$1 --bfile /rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/params/bed/chrX --load-model ${8}.fastGWA --geno 0.1 --out $9