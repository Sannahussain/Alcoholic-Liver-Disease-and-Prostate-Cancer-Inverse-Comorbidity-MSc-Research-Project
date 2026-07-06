Heritability:

# Estimate SNP heritability for PC
./ldsc.py \
  --h2 males_pc.sumstats.gz \
  --ref-ld-chr eur_w_ld_chr/ \
  --w-ld-chr eur_w_ld_chr/ \
  --out pc_h2

# Estimate SNP heritability for ALD
./ldsc.py \
  --h2 males_ald.sumstats.gz \
  --ref-ld-chr eur_w_ld_chr/ \
  --w-ld-chr eur_w_ld_chr/ \
  --out ald_h2
