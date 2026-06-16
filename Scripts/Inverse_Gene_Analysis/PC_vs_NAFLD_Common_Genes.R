#PC vs Non-alcoholic Fatty Liver Disease (NAFLD) common genes
#load Meta-analyses results for cross disease comparisons

# filter significant genes using FDR correction 
PCMetaA <- PCMetaA[PCMetaA$FDR < 0.05, ]

NAFLDMetaA <- NAFLDMetaA[NAFLDMetaA$FDR < 0.05, ]

#Find the overlap of common genes from Gene Symbol between PC and NAFLD
PCcommon_NAFLD <- intersect(PCMetaA$Gene_symbol , NAFLDMetaA$Gene_symbol )

#identify how many genes are shared between Meta-analyses 
length(PCcommon_NAFLD)

#subset the dataset to common genes, common genes from PC and then the common ones from NAFLD
PC_NAFLD <- PCMetaA[PCMetaA$Gene_symbol  %in% PCcommon_NAFLD, ]
NAFLD_common <- NAFLDMetaA[NAFLDMetaA$Gene_symbol  %in% PCcommon_NAFLD, ]

#order by Gene Symbol aligning rows
PC_NAFLD <- PC_NAFLD[order(PC_NAFLD$Gene_symbol ), ]
NAFLD_common <- NAFLD_common[order(NAFLD_common$Gene_symbol ), ]

#identify oppositely expressed genes
PCNAFLD_opposite_genes <- PC_NAFLD$Gene_symbol [
  (PC_NAFLD$Direction == "up" & NAFLD_common$Direction == "down") |
  (PC_NAFLD$Direction == "down" & NAFLD_common$Direction == "up")
]

#Identify how many genes are oppositely expressed
length(PCNAFLD_opposite_genes)

#results table
PCNAFLD_opposite_table <- data.frame(
   Gene_symbol  = PCNAFLD_opposite_genes,
  PC_Direction = PC_NAFLD$Direction[PC_NAFLD$Gene_symbol  %in% PCNAFLD_opposite_genes],
  NAFLD_Direction = NAFLD_common$Direction[NAFLD_common$Gene_symbol  %in% PCNAFLD_opposite_genes]
)

#view complete table
View(PCNAFLD_opposite_table)

#save file.csv
write.csv(PCNAFLD_opposite_table,
          "C:/Users/ahuss/Desktop/Research Project/Oppositely_Expressed_Genes/PCNAFLD_opposite_table.csv",
          row.names = FALSE
)
