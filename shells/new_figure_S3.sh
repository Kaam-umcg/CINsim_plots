#!/bin/bash

#SBATCH --job-name=new_Fig_S3
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16GB
#SBATCH --error=errors/new_Fig_S3.err
#SBATCH --output=outputs/new_Fig_S3.log

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
Rscript $TMPDIR/scripts/new_figure_S3.R

# copies all relevant results to SCRATCH
mkdir -p $SCRATCH/CINsim/new_Fig_S3
cp $TMPDIR $SCRATCH/CINsim/new_Fig_S3 -r