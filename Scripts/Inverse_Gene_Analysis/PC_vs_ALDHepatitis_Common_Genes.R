#PC vs ALDHepatitis common genes
#load Meta-analyses results for cross disease comparisons

# filter significant genes using FDR correction 
PCMetaA <- PCMetaA[PCMetaA$FDR < 0.05, ]

ALDhepMetaA <- ALDhepMetaA[ALDhepMetaA$FDR < 0.05, ]

#Find the overlap of common genes from Gene Symbol between PC and ALD hepatitis
PCcommon_ALDhep <- intersect(PCMetaA$Gene_symbol , ALDhepMetaA$Gene_symbol )

#identify how many genes are shared between Meta-analyses 
length(PCcommon_ALDhep)

#subset the dataset to common genes, common genes from PC and then the common ones from hep
PC_ALDhep <- PCMetaA[PCMetaA$Gene_symbol  %in% PCcommon_ALDhep, ]
ALDhep_common <- ALDhepMetaA[ALDhepMetaA$Gene_symbol  %in% PCcommon_ALDhep, ]

#order by Gene Symbol aligning rows
PC_ALDhep <- PC_ALDhep[order(PC_ALDhep$Gene_symbol ), ]
ALDhep_common <- ALDhep_common[order(ALDhep_common$Gene_symbol ), ]

#identify oppositely expressed genes
PChep_opposite_genes <- PC_ALDhep$Gene_symbol [
  (PC_ALDhep$Direction == "up" & ALDhep_common$Direction == "down") |
  (PC_ALDhep$Direction == "down" & ALDhep_common$Direction == "up")
]

#Identify how many genes are oppositely expressed
length(PChep_opposite_genes)

#results table
PChep_opposite_table <- data.frame(
   Gene_symbol  = PChep_opposite_genes,
  PC_Direction = PC_ALDhep$Direction[PC_ALDhep$Gene_symbol  %in% PChep_opposite_genes],
  ALDhep_Direction = ALDhep_common$Direction[ALDhep_common$Gene_symbol  %in% PChep_opposite_genes]
)

#View complete table
View(PChep_opposite_table)

#save file.csv
write.csv(PChep_opposite_table,
          "PChep_opposite_table.csv",
          row.names = FALSE
)
