#!/bin/bash
$1 --fastGWA-mlm --mbfile $2 --grm-sparse $3 --pheno $4 --mpheno $5 --qcovar $6 --covar $7 --out ${@:8}