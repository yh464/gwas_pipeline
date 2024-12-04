#!/bin/bash
source /home/yh464/.bashrc
conda activate gentoolspy

p=/rds/user/yh464/hpc-work/conda/polyfun/

python ${p}finemapper.py --method susie --n 54030 --sumstats $1 --chr $2 \
  --start $3 --end $4 --geno $5 --out $6 --max-num-causal 5 --allow-swapped-indel-alleles
