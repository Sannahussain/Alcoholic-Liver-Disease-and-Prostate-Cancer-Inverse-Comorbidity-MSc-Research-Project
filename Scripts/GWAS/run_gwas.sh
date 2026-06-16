#!/bin/bash

#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --array=1-22
#SBATCH --job-name=pc
#SBATCH --time=8:0:0

module add gcta

chr=$SLURM_ARRAY_TASK_ID
bgen=/path/to/ukb/Imputed/
grm=/path/to/grm/
covs=/path/to/covariates/

### males-only

## PC

gcta64 --fastGWA-mlm-binary \
--grm-sparse $grm/ea_sp_grm \
--pheno pc_pheno.txt \
--out pc/males_pc_${chr} \
--maf 0.001 \
--geno 0.02 \
--info 0.7 \
--keep males.txt \
--qcovar $covs/ageattend5pcs.txt \
--sample ukb22828_c22_b0_v3_s487164.sample \
--bgen $bgen/ukb_imp_chr${chr}_v3.bgen \
--thread-num ${SLURM_NTASKS}

## ALD

gcta64 --fastGWA-mlm-binary \
--grm-sparse $grm/ea_sp_grm \
--pheno ald_pheno.txt \
--out ald/males_ald_${chr} \
--maf 0.001 \
--geno 0.02 \
--info 0.7 \
--keep males.txt \
--qcovar $covs/ageattend5pcs.txt \
--sample ukb22828_c22_b0_v3_s487164.sample \
--bgen $bgen/ukb_imp_chr${chr}_v3.bgen \
--thread-num ${SLURM_NTASKS}

