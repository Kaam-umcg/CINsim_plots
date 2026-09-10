#!/bin/bash

#SBATCH --job-name=intro_plots
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --error=errors/intro_plots.err
#SBATCH --output=outputs/intro_plots.log

######
# MAKE SURE THIS IS THE CORRECT SCRIPT NAME
######
SCRIPT_NAME="intro_plots"
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