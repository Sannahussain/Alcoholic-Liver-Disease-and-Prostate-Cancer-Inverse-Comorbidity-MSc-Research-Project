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

## UK Biobank phenotype preprocessing
UK Biobank individual genotype and phenotype data were used for genetic epidemiology analyses, including GWAS, linkage disequilibrium score regression (LDSC) and local genetic correlation analysis (LAVA). Due to UK Biobank restriction, raw genotype amd participant level phenotype data cannot be included within this repository.
- Phenotype files for both ALD and PC were generated using R scripts wothin the HPC environment, included in the scripts repository. individuals were classified as cases or controls based on the predefined phenotype criteria. Binary phenotye files were generated, where affected individuals were assigned a value of 1 and unaffected were assigned a value of 0.
- Covariate files containing age and the first five genetic principal components (PCs) were generated for adjustment of population structure and potential confounding during GWAS analyses.

The following derived files were generated on the HPC environment but are not included in this repository due to data protection restrictions:
- `ald_pheno.txt` – ALD case-control phenotype file
- `pc_pheno.txt` – Prostate Cancer case-control phenotype file
- `ageattend5pcs.txt` – age and principal component covariate file
- `males.txt` – autosomal genotype sample file for male participants
- `males_X.txt` – X chromosome genotype sample file for male participants

# GWAS workflow
Genome wide association analyses were performed using male UK Biobank participants to reduce biological heterogeneity between sexes. Autosomal chromosomes (1–22) and chromosome X were analysed separately.

GWAS were performed using the fastGWA mixed linear model implemented in GCTA v1.94.1. Age and the first five genetic principal components were included as covariates to control for population stratification.

Analyses were performed using SLURM batch scripts on the HPC environment. Chromosome-specific summary statistics were merged and cleaned before downstream analyses.

# Linkage Disequilibrium Score Regression (LDSC)
- (LDSC analyses were performed using LDSC software (1.0.1) package in the HPC, cloned from the GitHub repository and executed using Python based scripts, including munge_sumstats.py for summary statistics pre-processing and ldsc.py for SNP heritability estimation. )
- 
LDSC was used to estimate SNP-based heritability and genome-wide genetic correlation between ALD and PC.

GWAS summary statistics were formatted using `munge_sumstats.py` before analysis. SNP heritability was estimated using univariate LDSC and genetic correlation was estimated using bivariate LDSC.

The European LD reference panel and HapMap3 SNP list (`w_hm3.snplist`) were used for LD score regression analyses.

Required LDSC resources:
- LDSC software package
- European LD score reference panel (`eur_w_ld_chr`)
- HapMap3 SNP list (`w_hm3.snplist`)

These reference files are not included in this repository and should be obtained from the LDSC repository 

# Local Genetic Correlation Analysis (LAVA)

LAVA was used to investigate local genetic correlation between ALD and PC across genomic regions.

The analysis used:
- GWAS summary statistics
- UK Biobank LD reference panel
- genomic locus definition file
- phenotype information file
- sample overlap matrix

The UK Biobank LD reference panel and locus definition file were generated/obtained separately and are not included due to file size and data restrictions.

