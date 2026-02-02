#!/bin/bash

#SBATCH --job-name=single_sweep_T_ALL
#SBATCH --time=5-00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=2GB
#SBATCH --error=errors/single_sweep_T_ALL.err
#SBATCH --output=outputs/single_sweep_T_ALL.log

######
# CHANGE THE NAME OF THE JOB TO THE SCRIPT NAME - OTHERWISE A LOT OF THINGS BREAK
######
START_TIME=$(date -d @"$SLURM_JOB_START_TIME")
END_TIME=$(date -d @"$SLURM_JOB_END_TIME")
echo "$SLURM_JOB_NAME" starting at "$START_TIME" with prospective end time of "$END_TIME"

# copy the script to the node
cp "../scripts/" "$TMPDIR" -r

# modules for running scripts
module purge
module load R/4.4.1-gfbf-2023b

# does the surv_FC also get over here?
echo "$1"
surv_FC="$1"
export surv_FC

# run your script
Rscript "$TMPDIR"/scripts/single_sweep_T302_p3.R

# copies all relevant results to SCRATCH
mkdir -p "$SCRATCH"/CINsim/single_sweep_T302_p3/"$surv_FC"
cp "$TMPDIR" "$SCRATCH"/CINsim/single_sweep_T302_p3/"$surv_FC" -r