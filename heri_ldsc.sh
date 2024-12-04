source /rds/user/yh464/hpc-work/ldsc-pyenv/bin/activate
ldsc_dir=/rds/user/yh464/hpc-work/ldsc/

cd /rds/user/yh464/rds-rb643-ukbiobank2/Data_Users/yh464/gwa
cd $1
if ! [ -d ../../gene_corr/ldsc_sumstats ]; then
  mkdir -p ../../gene_corr/ldsc_sumstats
fi

for x in *.fastGWA
do
  if ! [ -f ../../gene_corr/ldsc_sumstats/${x/fastGWA/sumstats} ]; then
    python ${ldsc_dir}munge_sumstats.py --sumstats $x --out ../../gene_corr/ldsc_sumstats/${x/.fastGWA/}
  fi
  if ! [ -f {x/fastGWA/h2.log} ]; then
    python ${ldsc_dir}ldsc.py --h2 ../../gene_corr/ldsc_sumstats/${x/fastGWA/sumstats} \
    --ref-ld-chr ${ldsc_dir}baseline/ --w-ld-chr ${ldsc_dir}baseline/ \
    --out ${x/fastGWA/h2}
  fi
done