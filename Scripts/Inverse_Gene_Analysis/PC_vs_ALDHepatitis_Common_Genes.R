#PC vs ALDHepatitis common genes
#load Meta-analyses results for cross disease comparisons

# filter significant genes using FDR correction 
PCMetaA <- PCMetaA[PCMetaA$FDR < 0.05, ]

ALDhepMetaA <- ALDhepMetaA[ALDhepMetaA$FDR < 0.05, ]

#common genes PC vs ALD hepatitis
PCcommon_ALDhep <- intersect(PCMetaA$Gene_symbol , ALDhepMetaA$Gene_symbol )
