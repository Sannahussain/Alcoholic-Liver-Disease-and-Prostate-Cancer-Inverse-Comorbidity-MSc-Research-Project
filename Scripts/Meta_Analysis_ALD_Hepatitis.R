#Meta-analysis of Differential Gene Expression in Alcoholic Hepatitis

#load datasets
ALDhepstudy1 <- read.csv("GSE143318_ALDhepatitis1signatureData.csv")
ALDhepstudy2 <- read.csv("GSE167308_ALDhepatitis5signatureData.csv")

#load required packages
library(metapod) 
library(dplyr) 

#find the overlap of Ensembl gene IDs across the datastets, these are the common genes between studies
ALDhepcommon_genes <- Reduce(intersect, list(ALDhepstudy1$Gene_symbol, ALDhepstudy2$Gene_symbol))

#identify how many genes are shared between all studies
length(ALDhepcommon_genes)

#remove genes not present in all studies so only common genes remain in each dataset
ALDhepstudy1 <- ALDhepstudy1[ALDhepstudy1$Gene_symbol %in% ALDhepcommon_genes, ] 
ALDhepstudy2 <- ALDhepstudy2[ALDhepstudy2$Gene_symbol %in% ALDhepcommon_genes, ] 


#reorder each dataset so they appear in the same order - so the datasets align correctly for meta-analysis
ALDhepstudy1 <- ALDhepstudy1[match(ALDhepcommon_genes, ALDhepstudy1$Gene_symbol), ]
ALDhepstudy2 <- ALDhepstudy2[match(ALDhepcommon_genes, ALDhepstudy2$Gene_symbol), ]


#store all the studies in a list for easier processing
ALDhepstudies <- list(ALDhepstudy1, ALDhepstudy2)

# Extract p-values and directions into separate lists
#convert these values for both lists to numerical values
ALDheppvals_list <- lapply(ALDhepstudies, function(x) as.numeric(as.character(x$PValue)))

# Extract numeric LogFC - indicates whether a gene is upregulated or downregulated
ALDheplogfc_list <- lapply(ALDhepstudies, function(x)
  as.numeric(as.character(x$Log_FoldChange)))


# Combine into long vectors
# combine p-values from all studies into one long vector
ALDhepall_pvals <- unlist(ALDheppvals_list)

#combine LogFC from all studies into one long vector
ALDhepall_logfc <- unlist(ALDheplogfc_list)

# Create grouping variable indicating which p-values belong to the same gene across different studies
ALDhepgene_group <- rep(ALDhepcommon_genes, times = length(ALDhepstudies))

#use the grouped slimes method to combine p-values for each gene, producing a single meta-analysis p-value per gene
ALDhepmeta_grouped <- groupedSimes(ALDhepall_pvals, ALDhepgene_group)

#determine the overall direction of change for each gene using the fold changes from the datasets
ALDhepgroup_dir <- summarizeGroupedDirection(
  effects = ALDhepall_logfc,
  group = ALDhepgene_group,
  influential = ALDhepmeta_grouped$influential
)

# Adjust for multiple testing - apply FDR corrrection
ALDhepfdr <- p.adjust(ALDhepmeta_grouped$p.value, method = "fdr")

# Final results table - combine gene IDs, meta-analysis p-values, adjusted p-values and direction of expresion change.
ALDhepresults <- data.frame(
  Gene_symbol = names(ALDhepmeta_grouped$p.value),
  Combined_P = ALDhepmeta_grouped$p.value,
  FDR = ALDhepfdr,
  Direction = ALDhepgroup_dir
)

#sort results by significance (lowest FDR first)
ALDhepresults <- ALDhepresults[order(ALDhepresults$FDR), ]

#view top results
head(ALDhepresults)

#open the complete table
View(ALDhepresults)

#check the total number of genes analysed
nrow(ALDhepresults)

#save the final meta-analysis results as a CSV file.
write.csv(ALDhepresults, 
          "C:/Users/ahuss/Desktop/Research Project/ALD Differential analysis/ALD_hepatitis_meta_results.csv",
          row.names = FALSE)

