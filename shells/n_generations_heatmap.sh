#!/bin/bash

#SBATCH --job-name=n_generations_heatmap
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16GB
#SBATCH --error=errors/n_generations_heatmap.err
#SBATCH --output=outputs/n_generations_heatmap.log

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

# the directory we'll be loading the data from
DATA_SOURCE="/scratch/p319788/CINsim/Figure_2_backup/tmp/results/T_ALL_params"
export DATA_SOURCE

# run your script
Rscript $TMPDIR/scripts/n_generations_heatmap.R

# copies all relevant results to SCRATCH
mkdir -p $SCRATCH/CINsim/n_generations_heatmap
cp $TMPDIR $SCRATCH/CINsim/n_generations_heatmap -r