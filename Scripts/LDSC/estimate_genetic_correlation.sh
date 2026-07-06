genetic correlation:

# Estimste genome-wide genetic correlation
./ldsc.py \
  --ref-ld-chr eur_w_ld_chr/ \
  --out pc_ald_rg \
  --rg males_pc.sumstats.gz,males_ald.sumstats.gz \
  --w-ld-chr eur_w_ld_chr/
