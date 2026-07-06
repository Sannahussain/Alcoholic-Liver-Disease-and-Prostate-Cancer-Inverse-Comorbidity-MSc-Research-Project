#LDSC (LD score regression) preprocessing:

#PC  preprocessing
munge_sumstats.py \
  --sumstats males_pc_all.fastGWA \
  --out PC \
  --N 210807

#ALD preprocessing
munge_sumstats.py \
  --sumstats males_ald_all.fastGWA \
  --out ald \
  --N 210807
