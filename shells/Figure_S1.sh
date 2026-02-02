#!/bin/bash

#SBATCH --job-name=Fig_S1
#SBATCH --time=0-12:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=16GB
#SBATCH --error=errors/Figure_S1.err
#SBATCH --output=outputs/Figure_S1.log

######
# CHANGE THE NAME OF THE JOB TO THE SCRIPT NAME - OTHERWISE A LOT OF THINGS BREAK
######
START_TIME=$(date -d @$SLURM_JOB_START_TIME)
END_TIME=$(date -d @$SLURM_JOB_END_TIME)
echo $SLURM_JOB_NAME starting at $START_TIME with prospective end time of $END_TIME

# copy the script to the node
cp "../scripts/" $TMPDIR -r

# modules for running scripts
module purge
module load R/4.4.1-gfbf-2023b

# run your script
Rscript $TMPDIR/scripts/Figure_S1.R

# copies all relevant results to SCRATCH
mkdir -p $SCRATCH/CINsim/Figure_S1
cp $TMPDIR $SCRATCH/CINsim/Figure_S1 -r