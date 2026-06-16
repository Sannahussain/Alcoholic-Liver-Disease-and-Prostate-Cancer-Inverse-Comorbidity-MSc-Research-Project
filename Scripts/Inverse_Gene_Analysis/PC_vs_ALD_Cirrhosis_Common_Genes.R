#PC vs ALDCirrhosis Common Genes
#load Meta-analysis PC results and Differential Gene Expression results for ALD Cirrhosis for cross disease comparisons

# filter significant genes using FDR correction 
PCMetaA <- PCMetaA[PCMetaA$FDR < 0.05, ]

ALDCirrhosisDE <- ALDCirrhosisDE[ALDCirrhosisDE$FDR_adjusted_P.value < 0.05, ]

#add Direction column to ALDCirrhosisDE dataset
ALDCirrhosisDE$Direction <- ifelse(ALDCirrhosisDE$Log_FoldChange > 0, "up", "down")

#Find the overlap of common genes from Gene Symbol between PC and ALD Cirrhosis
PCcommon_ALDCirr <- intersect(PCMetaA$Gene_symbol , ALDCirrhosisDE$Gene_symbol )

#identify how many genes are shared between both datasets
length(PCcommon_ALDCirr)

#subset the dataset to common genes, common genes from PC and then the common ones from Cirrhosis
PC_ALDCirr <- PCMetaA[PCMetaA$Gene_symbol %in% PCcommon_ALDCirr, ]
ALDCirr_common <- ALDCirrhosisDE[ALDCirrhosisDE$Gene_symbol %in% PCcommon_ALDCirr, ]

#order by Gene Symbol aligning rows
PC_ALDCirr <- PC_ALDCirr[order(PC_ALDCirr$Gene_symbol), ]
ALDCirr_common <- ALDCirr_common[order(ALDCirr_common$Gene_symbol), ]

#identify oppositely expressed genes
PCCirr_opposite_genes <- PC_ALDCirr$Gene_symbol[
  (PC_ALDCirr$Direction == "up" & ALDCirr_common$Direction == "down") |
  (PC_ALDCirr$Direction == "down" & ALDCirr_common$Direction == "up")
]

#Identify how many genes are oppositely expressed
length(PCCirr_opposite_genes)

#results table
PCCirr_opposite_table <- data.frame(
   Gene_symbol = PCCirr_opposite_genes,
  PC_Direction = PC_ALDCirr$Direction[PC_ALDCirr$Gene_symbol %in% PCCirr_opposite_genes],
  ALDCirr_Direction = ALDCirr_common$Direction[ALDCirr_common$Gene_symbol %in% PCCirr_opposite_genes]
)

#view complete table
View(PCCirr_opposite_table)

#save file.csv
write.csv(PCCirr_opposite_table,
          "PCCirr_opposite_table.csv",
          row.names = FALSE
)
