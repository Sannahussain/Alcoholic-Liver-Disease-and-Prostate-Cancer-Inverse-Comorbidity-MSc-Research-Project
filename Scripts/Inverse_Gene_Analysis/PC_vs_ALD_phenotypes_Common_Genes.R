#PC vs ALD phenotypes (Hepatitis and Cirrhosis) common genes

#Find the overlap of common genes from Gene Symbol between PC and ALD phenotypes

PC_Hep_Cirr_opposite <- Reduce(intersect, list(
  PCCirr_opposite_table$Gene_symbol,
  PChep_opposite_table$Gene_symbol
))

length(PC_Hep_Cirr_opposite)

#Subset each table to the common genes
Hep_PC_Cirr <- PChep_opposite_table[PChep_opposite_table$Gene_symbol %in% PC_Hep_Cirr_opposite, ]

Cirr_PC_Hep <- PCCirr_opposite_table[PCCirr_opposite_table$Gene_symbol %in% PC_Hep_Cirr_opposite, ]

#Order them to match
Hep_PC_Cirr <- Hep_PC_Cirr[order(Hep_PC_Cirr$Gene_symbol), ]
Cirr_PC_Hep <- Cirr_PC_Hep[order(Cirr_PC_Hep$Gene_symbol), ]


#Results table of common genes across all tables
PC_Hep_Cirr_opposite <- data.frame(
  Gene_symbol = Hep_PC_Cirr$Gene_symbol,
  PC = Hep_PC_Cirr$PC_Direction,
  Hepatitis = Hep_PC_Cirr$ALDhep_Direction,
  Cirrhosis = Cirr_PC_Hep$ALDCirr_Direction)

View(PC_Hep_Cirr_opposite)

#save file.csv
write.csv(PC_Hep_Cirr_opposite,
          "PC_Hep_Cirr_opposite.csv",
          row.names = FALSE
)

## reordering table to present - rank the common genes by overall significance and combined FDR values

#Rename column for clarity (ensure consistent naming of FDR column)
colnames(Cirr_FDR_lookup)[2] <- "FDR"

#Merge FDR results from different diseases
# Add prostate cancer FDR values to the combined dataset
PC_Hep_Cirr_opposite <- merge(PC_Hep_Cirr_opposite, PC_FDR_lookup,
                              by = "Gene_symbol", all.x = TRUE)

# Rename newly added column 
colnames(PC_Hep_Cirr_opposite)[ncol(PC_Hep_Cirr_opposite)] <- "PC_FDR"

# Add alcoholic hepatitis FDR values
PC_Hep_Cirr_opposite <- merge(PC_Hep_Cirr_opposite, Hep_FDR_lookup,
                              by = "Gene_symbol", all.x = TRUE)


colnames(PC_Hep_Cirr_opposite)[ncol(PC_Hep_Cirr_opposite)] <- "Hep_FDR"

# Add cirrhosis FDR values
PC_Hep_Cirr_opposite <- merge(PC_Hep_Cirr_opposite, Cirr_FDR_lookup,
                              by = "Gene_symbol", all.x = TRUE)

colnames(PC_Hep_Cirr_opposite)[ncol(PC_Hep_Cirr_opposite)] <- "Cirr_FDR"

#Create a combined significance score
# Convert FDR values into -log10 scale:
# Smaller FDR = more significant
# log transformation makes values easier to combine
# summing gives overall evidence across all diseases
PC_Hep_Cirr_opposite$combined_score <-
  -log10(PC_Hep_Cirr_opposite$PC_FDR) +
  -log10(PC_Hep_Cirr_opposite$Hep_FDR) +
  -log10(PC_Hep_Cirr_opposite$Cirr_FDR)

#Rank genes by overall significance
PC_Hep_Cirr_opposite <-
  PC_Hep_Cirr_opposite[
    order(-PC_Hep_Cirr_opposite$combined_score),
  ]

#view top 20 genes
head(PC_Hep_Cirr_opposite, 20)

#save file.csv
write.csv(PC_Hep_Cirr_opposite,
          "PC_Hep_Cirr_opposite.csv",
          row.names = FALSE
)
