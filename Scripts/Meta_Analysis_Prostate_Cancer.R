# Meta-analysis of Differential Gene Expression in Prostate Cancer

# Load datasets
PCstudy1 <- read.csv("GSE104131_PC1signatureData.csv") 
PCstudy2 <- read.csv("GSE22260_PC2signatureData.csv") 
PCstudy3 <- read.csv("GSE133626_PC3signatureData.csv")

#load required packages
library(metapod) 
library(dplyr) 


#find the overlap of Ensembl gene IDs across the datastets, these are the common genes between studies
PCcommon_genes <- Reduce(intersect, list(PCstudy1$Ensembl_ID, PCstudy2$Ensembl_ID, PCstudy3$Ensembl_ID)) 

#identify how many genes are shared between all studies
length(PCcommon_genes)

#remove genes not present in all studies so only common genes remain in each dataset
PCstudy1 <- PCstudy1[PCstudy1$Ensembl_ID %in% PCcommon_genes, ] 
PCstudy2 <- PCstudy2[PCstudy2$Ensembl_ID %in% PCcommon_genes, ] 
PCstudy3 <- PCstudy3[PCstudy3$Ensembl_ID %in% PCcommon_genes, ]


#reoder each dataset so they appear in the same order - so the datasets align correctly for meta-analysis
PCstudy1 <- PCstudy1[match(PCcommon_genes, PCstudy1$Ensembl_ID), ]
PCstudy2 <- PCstudy2[match(PCcommon_genes, PCstudy2$Ensembl_ID), ]
PCstudy3 <- PCstudy3[match(PCcommon_genes, PCstudy3$Ensembl_ID), ]


#store all the studies in a list for easier processing
PCstudies <- list(PCstudy1, PCstudy2, PCstudy3)

# Extract p-values and directions into separate lists
#convert these values for both lists to numerical values
PCpvals_list <- lapply(PCstudies, function(x) as.numeric(as.character(x$PValue)))

# Extract numeric LogFC - indicates whether a gene is upregulated or downregulated
PClogfc_list <- lapply(PCstudies, function(x)
  as.numeric(as.character(x$Log_FoldChange)))


# Combine into long vectors
#combine p-values from all studies into one long vector
PCall_pvals <- unlist(PCpvals_list)

#combine LogFC from all studies into one long vector
PCall_logfc <- unlist(PClogfc_list)

# Create grouping variable indicating which p-values belong to the same gene across different studies
PCgene_group <- rep(PCcommon_genes, times = length(PCstudies))

#use the grouped slimes method to combine p-values for each gene, producing a single meta-analysis p-value per gene
PCmeta_grouped <- groupedSimes(PCall_pvals, PCgene_group)

#determine the overall direction of change for each gene using the fold changes from the datasets
PCgroup_dir <- summarizeGroupedDirection(
  effects = PCall_logfc,
  group = PCgene_group,
  influential = PCmeta_grouped$influential
)

# Adjust for multiple testing - apply FDR correction
PCfdr <- p.adjust(PCmeta_grouped$p.value, method = "fdr")

# Final results table - combine gene IDs, meta-analysis p-values, adjusted p-values and direction of expresion change.
PCresults <- data.frame(
  Gene_symbol = names(PCmeta_grouped$p.value),
  Combined_P = PCmeta_grouped$p.value,
  FDR = PCfdr,
  Direction = PCgroup_dir
)

#sort results by significance (lowest FDR first)
PCresults <- PCresults[order(PCresults$FDR), ]

#view top results
head(PCresults)

#open the complete table
View(PCresults)

#check the total number of genes analysed
nrow(PCresults)

#save the final meta-analysis results as a CSV file.
write.csv(PCresults, "Prostate_cancer_meta_results.csv", row.names = FALSE)
