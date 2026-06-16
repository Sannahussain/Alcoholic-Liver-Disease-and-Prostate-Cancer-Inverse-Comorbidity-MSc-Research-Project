#create ALD phenotype file for GWAS

ald_pheno <- data.frame(
  FID = alcoset$id,
  IID = alcoset$id,
  ALD = alcoset$ald
)

#create ALD phenotype table
write.table(
  ald_pheno,
  "ald_pheno.txt",
  quote = FALSE,
  row.names = FALSE
)

#create PC phenotype file for GWAS
pc_pheno <- data.frame(
  FID = alcoset$id,
  IID = alcoset$id,
  PC = alcoset$pc
)

#create PC phenotype table
write.table(
  pc_pheno,
  "pc_pheno.txt",
  quote = FALSE,
  row.names = FALSE
)
