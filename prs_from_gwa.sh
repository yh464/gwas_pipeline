#!/bin/bash
#SBATCH -N 1 -n 1 -c 1 -t 12:00:00 -p cclake-himem
source /home/yh464/.bashrc
conda deactivate
conda activate gentoolspy

python ${1}/PRScs.py --ref_dir=${2} --bim_prefix=${3} --sst_file=${4} --n_gwas=${5} --out_dir=${6} --chrom=${7} --phi=${8} --seed 483647
