# SLURM submission script for LAVA local genetic correlation analysis
# This scrip submits the job array to the HPC to run the LAVA R script for each autosomal chromosome (1-22)
#==================================================================================================
#!/bin/bash

#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=pc_ald
#SBATCH --time=8:0:0

# Load the R module
module add r

# Run the LAVA analysis
Rscript lava_script.R "/path/to/ukb_reference/lava-ukb-v1.1" \
"/path/to/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile" \
"/path/to/input.info.txt" \
"/path/to/sample.overlap.txt" \
"PC;ALD" \
"pc_ald"
