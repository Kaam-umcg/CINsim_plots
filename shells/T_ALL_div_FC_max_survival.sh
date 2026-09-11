#!/bin/bash

#SBATCH --job-name=T_ALL_div_FC_max_survival
#SBATCH --time=5-00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4GB
#SBATCH --error=errors/T_ALL_div_FC_max_survival.err
#SBATCH --output=outputs/T_ALL_div_FC_max_survival.log

######
# MAKE SURE THIS IS THE CORRECT SCRIPT NAME
######
SCRIPT_NAME="T_ALL_div_FC_max_survival"
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

# run your script
Rscript "$TMPDIR/$SCRIPT_NAME.R"

# copies all relevant results to SCRATCH
mkdir -p "$SCRATCH/CINsim/$SCRIPT_NAME"
cp "$TMPDIR"/* "$SCRATCH"/CINsim/$SCRIPT_NAME -r