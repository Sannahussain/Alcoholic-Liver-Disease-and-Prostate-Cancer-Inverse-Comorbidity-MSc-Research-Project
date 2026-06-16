#PC vs all liver disease phenotypes common genes

#Find the overlap of common genes from Gene Symbol between PC and all liver disease phenotypes
common_all_opposite <- Reduce(intersect, list(
  PCCirr_opposite_table$Gene_symbol,
  PCNASH_opposite_table$Gene_symbol,
  PCNAFLD_opposite_table$Gene_symbol,
  PChep_opposite_table$Gene_symbol
))

#identify how many genes are shared between Datasets
length(common_all_opposite)

#subset the datasets to common genes, common genes from PC and then the common ones from eahc liver disease phenotype
Hep_sub <- PChep_opposite_table[PChep_opposite_table$ Gene_symbol %in% common_all_opposite, ]

Cirr_sub <- PCCirr_opposite_table[PCCirr_opposite_table$Gene_symbol %in% common_all_opposite, ]

NAFLD_sub <- PCNAFLD_opposite_table[PCNAFLD_opposite_table$Gene_symbol %in% common_all_opposite, ]

NASH_sub <- PCNASH_opposite_table[PCNASH_opposite_table$Gene_symbol %in% common_all_opposite, ]

#order by Gene Symbol aligning rows
Hep_sub <- Hep_sub[order(Hep_sub$Gene_symbol), ]
Cirr_sub <- Cirr_sub[order(Cirr_sub$Gene_symbol), ]
NAFLD_sub <- NAFLD_sub[order(NAFLD_sub$Gene_symbol), ]
NASH_sub <- NASH_sub[order(NASH_sub$Gene_symbol), ]

#Results table of common genes across all tables
common_all_opposite_table <- data.frame(
  Gene_symbol = Hep_sub$Gene_symbol,
  PC = Hep_sub$PC_Direction,
  Hepatitis = Hep_sub$ALDhep_Direction,
  Cirrhosis = Cirr_sub$ALDCirr_Direction,
  NAFLD = NAFLD_sub$NAFLD_Direction,
  NASH = NASH_sub$NASH_Direction
)

#view complete table
View(common_all_opposite_table)

#save file.csv
write.csv(common_all_opposite_table,
          "common_all_opposite_table.csv",
          row.names = FALSE
)
