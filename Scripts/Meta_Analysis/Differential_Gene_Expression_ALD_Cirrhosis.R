# Differential Gene Expression in Alcoholic Cirrhosis Table

#Load study
Cirrhosisstudy1 <- read.csv("GSE167308_ALDcirrhosis5signatureData.csv" )

#identify how many genes present in the study
length(Cirrhosisstudy1)

#view genes
head(Cirrhosisstudy1)

#open complete dataset
View(Cirrhosisstudy1)

#save the Differential Gene expression table
write.csv(Cirrhosisstudy1, 
          "ALD Differential analysis/ALD_Cirrhosis_DE_results.csv",
          row.names = FALSE)
