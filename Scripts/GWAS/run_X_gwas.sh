#!/bin/bash

#SBATCH --mem=30G
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --job-name=pc
#SBATCH --time=4:0:0

module add gcta

module add gcta

chr=$SLURM_ARRAY_TASK_ID
bgen=/path/to/ukb/Imputed/
grm=/path/to//grm/
covs=/path/to/covariates/

### males-only

## PC
# model
gcta64 --fastGWA-mlm-binary \
--grm-sparse $grm/ea_sp_grm \
--pheno pc_pheno.txt \
--out pcX/males_model_X \
--qcovar $covs/ageattend5pcs.txt \
--sample ukb22828_c22_b0_v3_s487164.sample \
--keep males_X.txt \
--bgen $bgen/ukb_imp_chr22_v3.bgen \
--model-only \
--thread-num ${SLURM_NTASKS}

# X
gcta64 --bgen $bgen/ukb_imp_chrX_v3.bgen \
--sample ukb22828_cX_b0_v3_s486511.sample \
--load-model pcX/males_model_X.fastGWA \
--keep males_X.txt \
--maf 0.001 --geno 0.02 --info 0.7 \
--out pcX/males_pc_23 \
--thread-num ${SLURM_NTASKS}

####### Same for ALD
gcta64 --fastGWA-mlm-binary \
--grm-sparse $grm/ea_sp_grm \
--pheno ald_pheno.txt \
--out aldX/males_model_X \
--qcovar $covs/ageattend5pcs.txt \
--sample ukb22828_c22_b0_v3_s487164.sample \
--keep males_X.txt \
--bgen $bgen/ukb_imp_chr22_v3.bgen \
--model-only \
--thread-num ${SLURM_NTASKS}

gcta64 --bgen $bgen/ukb_imp_chrX_v3.bgen \
--sample ukb22828_cX_b0_v3_s486511.sample \
--load-model aldX/males_model_X.fastGWA \
--keep males_X.txt \
--maf 0.001 --geno 0.02 --info 0.7 \
--out aldX/males_ald_23 \
--thread-num ${SLURM_NTASKS}

