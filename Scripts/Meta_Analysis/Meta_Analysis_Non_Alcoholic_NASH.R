#Meta-analysis of Differential Gene Expression in Non-Alcoholic Steatohepatitis (NASH)

#load datasets
NASHstudy1 <- read.csv("GSE126848_NON-ALD2_NASHvsHEALTHY_signatureData.csv") 
NASHstudy2 <- read.csv("GSE260666_NON-ALD_NASH4signatureData.csv") 

#load rewuired packages
library(metapod) 
library(dplyr) 


#find the overlap of Ensembl gene symbols across the datastets, these are the common genes between studies
NASHcommon_genes <- Reduce(intersect, list(NASHstudy1$Gene_symbol, NASHstudy2$Gene_symbol)) 

#identify how many genes are shared between all studies
length(NASHcommon_genes)

#remove genes not present in all studies so only common genes remain in each dataset
NASHstudy1 <- NASHstudy1[NASHstudy1$Gene_symbol %in% NASHcommon_genes, ] 
NASHstudy2 <- NASHstudy2[NASHstudy2$Gene_symbol %in% NASHcommon_genes, ] 

#reorder each dataset so they appear in the same order - so the datasets align correctly for meta-analysis
NASHstudy1 <- NASHstudy1[match(NASHcommon_genes, NASHstudy1$Gene_symbol), ]
NASHstudy2 <- NASHstudy2[match(NASHcommon_genes, NASHstudy2$Gene_symbol), ]

#store all the studies in a list for easier processing
NASHstudies <- list(NASHstudy1, NASHstudy2)


# Extract p-values and directions into separate lists
#convert these values for both lists to numerical values
NASHpvals_list <- lapply(NASHstudies, function(x) as.numeric(as.character(x$PValue)))

# Extract numeric LogFC - indicates whether a gene is upregulated or downregulated
NASHlogfc_list <- lapply(NASHstudies, function(x)
  as.numeric(as.character(x$Log_FoldChange)))


# Combine into long vectors
# combine p-values from all studies into one long vector
NASHall_pvals <- unlist(NASHpvals_list)

#combine LogFC from all studies into one long vector
NASHall_logfc <- unlist(NASHlogfc_list)


# Create grouping variable indicating which p-values belong to the same gene across different studies
NASHgene_group <- rep(NASHcommon_genes, times = length(NASHstudies))

#use the grouped slimes method to combine p-values for each gene, producing a single meta-analysis p-value per gene
NASHmeta_grouped <- groupedSimes(NASHall_pvals, NASHgene_group)

#determine the overall direction of change for each gene using the fold changes from the datasets
NASHgroup_dir <- summarizeGroupedDirection(
  effects = NASHall_logfc,
  group = NASHgene_group,
  influential = NASHmeta_grouped$influential
)

# Adjust for multiple testing - apply FDR correction
NASHfdr <- p.adjust(NASHmeta_grouped$p.value, method = "fdr")


# Final results table - combine gene symbols, meta-analysis p-values, adjusted p-values and direction of expresion change.
NASHresults <- data.frame(
  Gene_symbol = names(NASHmeta_grouped$p.value),
  Combined_P = NASHmeta_grouped$p.value,
  FDR = NASHfdr,
  Direction = NASHgroup_dir
)

#sort results by significance (lowest FDR first)
NASHresults <- NASHresults[order(NASHresults$FDR), ]

#view top results
head(NASHresults)

#open the complete table
View(NASHresults)

#check the total number of genes analysed
nrow(NASHresults)

#save the final meta-analysis results as a CSV file.
write.csv(NASHresults, 
          "NON-ALD_NASH_meta_results.csv",
          row.names = FALSE)
