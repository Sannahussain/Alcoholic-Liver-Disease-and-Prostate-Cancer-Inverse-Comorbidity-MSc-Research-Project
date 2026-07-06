# Alcoholic-Liver-Disease-and-Prostate-Cancer-Inverse-Comorbidity-MSc-Research-Project
PROJECT TITLE:
The Molecular Genetic Basis of Inverse Comorbidity of Alcoholic Liver Disease (ALD) and Prostate Cancer (PC).

# Overview:
This repository contains scripts used to perform the transcriptomic and genetic epidemiology analyses between PC and ALD to investigating their genetic relationship, including genome wode and post GWAS analyses.

# The workflow includes:
- Meta analysis of liver disease and prostate cancer phenotypes
- Inverse gene level comparisons between liver disease and prostae cancer phenotypes
- Genome Wide Association Studies (GWAS)
- Manhattan plot visualisation
- SNP heritibility estimation (LDSC)
- Genetic correlation analysis (LDSC)
- Local genetic correlation analysis (LAVA)

#  Software Requirements:
- R
- PLINK
- LDSC
- LAVA
- Python (for LDSC preprocessing)
- HPC environment

# Data Sources:
Differemtial gene expression and meta analyses were performed using publicly available RNA sequencing datasets from the GREIN database.
# Prostate Cancer datasets
  - Dataset GSE104131 - RNA sequencing of prostate cancer and matched normal tissue from African American and European American individuals (16 tumour, 16 control samples)
  - Dataset GSE22260 - Comparative transcriptomic analysis of prostate cancer and matched normal tissue (20 tumour, 10 control samples)
  - Dataset GSE133626	- Patient matched prostate cancer and normal tissue RNA sequencing (30 tumour, 30 control samples)

# Liver Disease datasets
  - Dataset GSE142218 - Alcoholic hepatitis (5 cases, 5 controls)
  - Dataset GSE167308 - Alcoholic cirrhosis and alcoholic hepatitis (7 cirrhosis, 7 hepatitis, 5 controls)
  - Dataset GSE126848 - Non-alcoholic fatty liver disease (NAFLD) and non-alcoholic steatohepatitis (NASH) (15 NAFLD, 16 NASH, 14 controls) 
  - Dataset GSE260666 - NAFLD and NASH (6 NAFLD, 4 NASH, 6 controls)

All RNA-seq datasets were used to generate separate signature IDs for each phenotype for affected vs control samples before being used in downstream meta analysis and inverse gene level comparison workflows.

 # UK Biobank
UK Biobank data, which is not publicly available due to access restriction, was used for genetic epidemiology analyses including GWAS, LDSC and LAVA analyses.
- Summary statistics and intermediate files are generated on HPC and are not included in this repository. 




