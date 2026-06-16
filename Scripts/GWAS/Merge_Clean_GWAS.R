#Merge and Clean GWAS

# remove columns not required

#load required library
library(data.table)

#PC cleaning
files <- list.files("/path/to/pc/.fastGWA/",
                    pattern = "males_pc_.*\\.fastGWA$",
                    full.names = TRUE)

for (f in files) {

  df <- fread(f)

  df_clean <- df[, .(
    CHR, SNP, POS, A1, A2, N, AF1, BETA, P
  )]

  out <- gsub("males_", "clean_males_", f)

  fwrite(df_clean, out, sep = "\t")
}

#PC merging cleaned .fastGWA files

files <- list.files("/path/to/clean/pc/.fastGWA",
                    pattern = "clean_males_pc_.*\\.fastGWA$",
                    full.names = TRUE)

pc_list <- lapply(files, fread)

males_pc_all <- rbindlist(pc_list)

fwrite(males_pc_all,
       "/path/to/merged/males_pc_all.fastGWA",
       sep = "\t")


#ALD cleaning

files <- list.files("/path/to/ald/.fastGWA/",
                    pattern = "males_ald_.*\\.fastGWA$",
                    full.names = TRUE)

for (f in files) {

  df <- fread(f)

  df_clean <- df[, .(
    CHR, SNP, POS, A1, A2, N, AF1, BETA, P
  )]

  out <- gsub("males_", "clean_males_", f)

  fwrite(df_clean, out, sep = "\t")
}

#ALD merging cleaned .fastGWA files

files <- list.files("/path/to/ald/.fastGWA",
                    pattern = "clean_males_ald_.*\\.fastGWA$",
                    full.names = TRUE)

ald_list <- lapply(files, fread)

males_ald_all <- rbindlist(ald_list)

fwrite(males_ald_all,
       "/path/to/merged/males_ald_all.fastGWA",
       sep = "\t")
