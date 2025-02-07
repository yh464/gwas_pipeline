#!/bin/bash
source /home/yh464/.bashrc
conda deactivate
conda activate $1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
Rscript ${@:2}