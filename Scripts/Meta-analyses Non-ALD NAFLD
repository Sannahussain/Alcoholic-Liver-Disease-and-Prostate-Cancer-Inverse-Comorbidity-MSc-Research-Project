#Meta-analysis of Differential Gene Expression in Non-Alcoholic Fatty Liver Disease (NAFLD)

#load datasets
NAFLDstudy1 <- read.csv("GSE126848_NON-ALD2_NAFLDvsHEALTHY_signatureData.csv") 
NAFLDstudy2 <- read.csv("GSE260666_NON-ALD_NAFLD4signatureData.csv") 

#load required packages
library(metapod) 
library(dplyr) 

#find the overlap of Ensembl gene IDs across the datastets, these are the common genes between studies
NAFLDcommon_genes <- Reduce(intersect, list(NAFLDstudy1$Gene_symbol, NAFLDstudy2$Gene_symbol)) 

#identify how many genes are shared between all studies
length(NAFLDcommon_genes)

#remove genes not present in all studies so only common genes remain in each dataset
NAFLDstudy1 <- NAFLDstudy1[NAFLDstudy1$Gene_symbol %in% NAFLDcommon_genes, ] 
NAFLDstudy2 <- NAFLDstudy2[NAFLDstudy2$Gene_symbol %in% NAFLDcommon_genes, ] 

#reorder each dataset so they appear in the same order - so the datasets align correctly for meta-analysis
NAFLDstudy1 <- NAFLDstudy1[match(NAFLDcommon_genes, NAFLDstudy1$Gene_symbol), ]
NAFLDstudy2 <- NAFLDstudy2[match(NAFLDcommon_genes, NAFLDstudy2$Gene_symbol), ]

#store all the studies in a list for easier processing
NAFLDstudies <- list(NAFLDstudy1, NAFLDstudy2)


# Extract p-values and directions into separate lists
#convert these values for both lists to numerical values
NAFLDpvals_list <- lapply(NAFLDstudies, function(x) as.numeric(as.character(x$PValue)))

# Extract numeric LogFC - indicates whether a gene is upregulated or downregulated
NAFLDlogfc_list <- lapply(NAFLDstudies, function(x)
  as.numeric(as.character(x$Log_FoldChange)))


# Combine into long vectors
# combine p-values from all studies into one long vector
NAFLDall_pvals <- unlist(NAFLDpvals_list)

#combine LogFC from all studies into one long vector
NAFLDall_logfc <- unlist(NAFLDlogfc_list)

# Create grouping variable indicating which p-values belong to the same gene across different studies
NAFLDgene_group <- rep(NAFLDcommon_genes, times = length(NAFLDstudies))

#use the grouped slimes method to combine p-values for each gene, producing a single meta-analysis p-value per gene
NAFLDmeta_grouped <- groupedSimes(NAFLDall_pvals, NAFLDgene_group)

#determine the overall direction of change for each gene using the fold changes from the datasets
NAFLDgroup_dir <- summarizeGroupedDirection(
  effects = NAFLDall_logfc,
  group = NAFLDgene_group,
  influential = NAFLDmeta_grouped$influential
)

# Adjust for multiple testing - apply FDR correction
NAFLDfdr <- p.adjust(NAFLDmeta_grouped$p.value, method = "fdr")


# Final results tablee - combine gene IDs, meta-analysis p-values, adjusted p-values and direction of expresion change.
NAFLDresults <- data.frame(
  Gene_symbol = names(NAFLDmeta_grouped$p.value),
  Combined_P = NAFLDmeta_grouped$p.value,
  FDR = NAFLDfdr,
  Direction = NAFLDgroup_dir
)

#sort results by significance (lowest FDR first)
NAFLDresults <- NAFLDresults[order(NAFLDresults$FDR), ]

#view top results
head(NAFLDresults)

#open the complete table
View(NAFLDresults)

#check the total number of genes analysed
nrow(NAFLDresults)

#save the final meta-analysis results as a CSV file.
write.csv(NAFLDresults, 
          "C:/Users/ahuss/Desktop/Research Project/ALD Differential analysis/NON-ALD_NAFLD_meta_results.csv",
          row.names = FALSE)
