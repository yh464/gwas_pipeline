#!/bin/bash
$1 --reml --reml-est-fix --reml-pred-rand --grm-cutoff 0.05 --grm $2 --pheno $3 --mpheno $4 --covar $5 --qcovar $6 --out ${@:7}
