#!/bin/bash

#SBATCH --job-name=organoid_CnFS_barplot
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8GB
#SBATCH --error=errors/organoid_CnFS_barplot.err
#SBATCH --output=outputs/organoid_CnFS_barplot.log

######
# MAKE SURE THIS IS THE CORRECT SCRIPT NAME
######
SCRIPT_NAME="organoid_CnFS_barplot"
######
# MAKE SURE THIS IS THE CORRECT SCRIPT NAME
######

START_TIME=$(date -d @"$SLURM_JOB_START_TIME")
END_TIME=$(date -d @"$SLURM_JOB_END_TIME")
echo "$SLURM_JOB_NAME.R" starting at "$START_TIME" with prospective end time of "$END_TIME"

# copy the script to the node
cp "../scripts/$SCRIPT_NAME.R" "$TMPDIR"

# modules for running scripts
module purge
module load R/4.4.1-gfbf-2023b

#TODO update these
export "BASE_PATH"="/scratch/p319788/CINsim"
export "DIPLOID_PATH"="Figure_S4a_diploid/tmp/results"
export "WGD_PATH"="Figure_S4a_WGD/tmp/results"
export "OUTPUT_PATH"="/scratch/p319788/CINsim/Figure_S4a_plots"

# run your script
Rscript "$TMPDIR/$SCRIPT_NAME.R"

# copies all relevant results to SCRATCH
mkdir -p "$SCRATCH/CINsim/$SCRIPT_NAME"
cp "$TMPDIR"/* "$SCRATCH"/CINsim/$SCRIPT_NAME -r
