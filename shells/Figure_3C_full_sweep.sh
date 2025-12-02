#!/bin/bash

#SBATCH --job-name=full_sweep_T_ALL
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --error=errors/full_sweep_T_ALL.err
#SBATCH --output=outputs/full_sweep_T_ALL.log

######
# CHANGE THE NAME OF THE JOB TO THE SCRIPT NAME - OTHERWISE A LOT OF THINGS BREAK
######
START_TIME=$(date -d @"$SLURM_JOB_START_TIME")
END_TIME=$(date -d @"$SLURM_JOB_END_TIME")
echo "$SLURM_JOB_NAME" starting at "$START_TIME" with prospective end time of "$END_TIME"

for surv_FC in 1.11 1.15 1.20 1.25 1.30 1.36 1.43 1.50 1.58 1.67 1.76 1.88 2.00
do
    echo "Starting script for survival_FC of $surv_FC"  

    # passes the surv_FC in the sbatch call
    sbatch --output="outputs/T_ALL_single_sweep_$surv_FC.log" --error="errors/T_ALL_single_sweep_$surv_FC.err" single_sweep_T_ALL.sh "$surv_FC"

    # sleep to cancel in case something goes wrong
    sleep 120
done
