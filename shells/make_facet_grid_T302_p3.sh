#!/bin/bash

#SBATCH --job-name=make_facet_grid_T302_p3
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8GB
#SBATCH --error=errors/make_facet_grid_T302_p3.err
#SBATCH --output=outputs/make_facet_grid_T302_p3.log

######
# CHANGE THE NAME OF THE JOB TO THE SCRIPT NAME - OTHERWISE A LOT OF THINGS BREAK
######
START_TIME=$(date -d @"$SLURM_JOB_START_TIME")
END_TIME=$(date -d @"$SLURM_JOB_END_TIME")
echo "$SLURM_JOB_NAME" starting at "$START_TIME" with prospective end time of $END_TIME

# copy the script to the node
cp "../scripts/make_facet_grid_T302_p3.R" "$TMPDIR"

# modules for running scripts
module purge
module load R/4.4.1-gfbf-2023b

# run your script
Rscript "$TMPDIR"/make_facet_grid_T302_p3.R

# copies all relevant results to SCRATCH
mkdir -p "$SCRATCH"/CINsim/make_facet_grid_T302_p3
cp "$TMPDIR" "$SCRATCH"/CINsim/make_facet_grid_T302_p3 -r
