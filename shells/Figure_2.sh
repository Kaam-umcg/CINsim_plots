#!/bin/bash

#SBATCH --job-name=Figure_2
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2GB
#SBATCH --error=errors/Figure_2.err
#SBATCH --output=outputs/Figure_2.log

######
# CHANGE THE NAME OF THE JOB TO THE SCRIPT NAME - OTHERWISE A LOT OF THINGS BREAK
######
START_TIME=$(date -d @$SLURM_JOB_START_TIME)
END_TIME=$(date -d @$SLURM_JOB_END_TIME)
echo $SLURM_JOB_NAME starting at $START_TIME with prospective end time of $END_TIME

# copy the script to the node
cp "../scripts/Figure_2.R" $TMPDIR 

# modules for running scripts
module purge
module load R/4.4.1-gfbf-2023b

# run your script
Rscript $TMPDIR/Figure_2.R

# copies all relevant results to SCRATCH
mkdir -p $SCRATCH/CINsim/Figure_2 
cp $TMPDIR $SCRATCH/CINsim/Figure_2 -r