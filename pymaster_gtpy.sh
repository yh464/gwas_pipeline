#!/bin/bash
source /home/yh464/.bashrc
module load gcc/11
# conda activate wd
source /rds/user/yh464/hpc-work/yh464-pyenv/bin/activate
conda activate gentoolspy
python /rds/project/rb643-1/rds-rb643-ukbiobank2/Data_Users/yh464/scripts/$@
