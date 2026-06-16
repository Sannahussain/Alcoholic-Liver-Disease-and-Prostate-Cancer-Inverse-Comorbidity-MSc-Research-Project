#PC vs Non-alcoholic Steatohepatitis (NASH) common genes
#load Meta-analyses results for cross disease comparisons

# filter significant genes using FDR correction 
PCMetaA <- PCMetaA[PCMetaA$FDR < 0.05, ]

NASHMetaA <- NASHMetaA[NASHMetaA$FDR < 0.05, ]

#Find the overlap of common genes from Gene Symbol between PC and NAFLD
PCcommon_NASH <- intersect(PCMetaA$Gene_symbol , NASHMetaA$Gene_symbol )

#identify how many genes are shared between Meta-analyses 
length(PCcommon_NASH)

#subset the dataset to common genes, common genes from PC and then the common ones from NASH
PC_NASH <- PCMetaA[PCMetaA$Gene_symbol  %in% PCcommon_NASH, ]
NASH_common <- NASHMetaA[NASHMetaA$Gene_symbol  %in% PCcommon_NASH, ]


#order by Gene Symbol aligning rows
PC_NASH <- PC_NASH[order(PC_NASH$Gene_symbol ), ]
NASH_common <- NASH_common[order(NASH_common$Gene_symbol ), ]

#identify oppositely expressed genes
PCNASH_opposite_genes <- PC_NASH$Gene_symbol [
  (PC_NASH$Direction == "up" & NASH_common$Direction == "down") |
  (PC_NASH$Direction == "down" & NASH_common$Direction == "up")
]

#Identify how many genes are oppositely expressed
length(PCNASH_opposite_genes)

#results table
PCNASH_opposite_table <- data.frame(
   Gene_symbol  = PCNASH_opposite_genes,
  PC_Direction = PC_NASH$Direction[PC_NASH$Gene_symbol  %in% PCNASH_opposite_genes],
  NASH_Direction = NASH_common$Direction[NASH_common$Gene_symbol  %in% PCNASH_opposite_genes]
)

#view complete table
View(PCNASH_opposite_table)

#save file.csv
write.csv(PCNASH_opposite_table,
          "C:/Users/ahuss/Desktop/Research Project/Oppositely_Expressed_Genes/PCNASH_opposite_table.csv",
          row.names = FALSE
)
